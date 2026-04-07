#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

#include "../src/payoffs/MultiAssetPayoff.hpp"
#include "../src/pricing/MultiAssetAnalytical.hpp"
#include "../src/pricing/MultiAssetPricer.hpp"
#include "../src/simulators/MultiAssetQMCSimulator.hpp"
#include "../src/stochasticProcess/CorrelatedGBM.hpp"

namespace {

const std::string OUTPUT_DIR = "tests/output";
const std::string CSV_FILE = OUTPUT_DIR + "/multiasset_ordering_convergence.csv";
const std::string PLOT_SCRIPT = OUTPUT_DIR + "/plot_multiasset_ordering_convergence.py";

const std::vector<int> PATH_COUNTS = {256, 512, 1024, 2048, 4096};
const int MC_OUTER_RUNS = 4;

struct SeriesPoint {
    std::string label;
    int paths;
    double error;
    double errorStdDev;
    double slope;
};

double mean(const std::vector<double>& xs) {
    return xs.empty() ? 0.0 : std::accumulate(xs.begin(), xs.end(), 0.0) / xs.size();
}

double stdDev(const std::vector<double>& xs) {
    if (xs.size() < 2)
        return 0.0;
    const double mu = mean(xs);
    double acc = 0.0;
    for (double x : xs) {
        const double d = x - mu;
        acc += d * d;
    }
    return std::sqrt(acc / (xs.size() - 1));
}

double slopeFromLogLog(const std::vector<int>& x, const std::vector<double>& y) {
    std::vector<double> xs;
    std::vector<double> ys;
    for (std::size_t i = 0; i < x.size(); ++i) {
        if (x[i] > 0 && y[i] > 0.0 && std::isfinite(y[i])) {
            xs.push_back(std::log(static_cast<double>(x[i])));
            ys.push_back(std::log(y[i]));
        }
    }
    if (xs.size() < 2)
        return std::numeric_limits<double>::quiet_NaN();

    const double mx = mean(xs);
    const double my = mean(ys);
    double num = 0.0;
    double den = 0.0;
    for (std::size_t i = 0; i < xs.size(); ++i) {
        num += (xs[i] - mx) * (ys[i] - my);
        den += (xs[i] - mx) * (xs[i] - mx);
    }
    return den > 0.0 ? num / den : std::numeric_limits<double>::quiet_NaN();
}

void ensureOutputDir() {
#ifdef _WIN32
    _mkdir(OUTPUT_DIR.c_str());
#else
    mkdir(OUTPUT_DIR.c_str(), 0755);
#endif
}

double priceFromSimulation(
    const MultiAssetQMCSimulator::SimulationResult& sim,
    const MultiAssetPayoff& payoff,
    double r,
    double T)
{
    return MultiAssetPricer::price(sim.terminalPrices, payoff, r, T).first;
}

SeriesPoint runMcPoint(
    const CorrelatedGBM& process,
    const std::vector<double>& s0,
    const MultiAssetPayoff& payoff,
    double r,
    double T,
    double dt,
    int paths,
    double exactPrice)
{
    std::vector<double> prices;
    prices.reserve(MC_OUTER_RUNS);
    for (int run = 0; run < MC_OUTER_RUNS; ++run) {
        const auto sim = MultiAssetQMCSimulator::simulate(
            process, s0, T, dt, paths,
            MultiAssetQMCSimulator::MC,
            "docs/new-joe-kuo-6.21201.txt",
            true, false, 8,
            static_cast<std::uint64_t>(42 + run * 101),
            MultiAssetQMCSimulator::Cholesky,
            MultiAssetQMCSimulator::AssetFirst);
        prices.push_back(priceFromSimulation(sim, payoff, r, T));
    }

    std::vector<double> absErr(prices.size(), 0.0);
    double sq = 0.0;
    for (std::size_t i = 0; i < prices.size(); ++i) {
        absErr[i] = std::abs(prices[i] - exactPrice);
        const double d = prices[i] - exactPrice;
        sq += d * d;
    }

    SeriesPoint point;
    point.label = "MC";
    point.paths = paths;
    point.error = std::sqrt(sq / prices.size());
    point.errorStdDev = stdDev(absErr);
    point.slope = std::numeric_limits<double>::quiet_NaN();
    return point;
}

SeriesPoint runQmcPoint(
    const std::string& label,
    const CorrelatedGBM& process,
    const std::vector<double>& s0,
    const MultiAssetPayoff& payoff,
    double r,
    double T,
    double dt,
    int paths,
    double exactPrice,
    MultiAssetQMCSimulator::PathConstruction construction,
    MultiAssetQMCSimulator::DimensionOrdering ordering)
{
    const auto sim = MultiAssetQMCSimulator::simulate(
        process, s0, T, dt, paths,
        MultiAssetQMCSimulator::QMC,
        "docs/new-joe-kuo-6.21201.txt",
        true, false, 8, 42,
        construction, ordering, &payoff, 24);
    const double price = priceFromSimulation(sim, payoff, r, T);

    SeriesPoint point;
    point.label = label;
    point.paths = paths;
    point.error = std::abs(price - exactPrice);
    point.errorStdDev = 0.0;
    point.slope = std::numeric_limits<double>::quiet_NaN();
    return point;
}

void writePlotScript() {
    std::ofstream py(PLOT_SCRIPT.c_str());
    py << "import csv\n";
    py << "from collections import defaultdict\n";
    py << "from pathlib import Path\n\n";
    py << "csv_path = Path(r'" << CSV_FILE << "')\n";
    py << R"PY(
rows = list(csv.DictReader(csv_path.open(newline='')))
out_svg = csv_path.with_suffix('.svg')
colors = {
    'MC': '#222222',
    'QMC asset-first': '#1f77b4',
    'QMC time-first': '#ff7f0e',
    'QMC PCA+bridge': '#2ca02c',
    'QMC adaptive': '#d62728',
}
series = defaultdict(lambda: {'x': [], 'y': [], 'slope': None})
for row in rows:
    series[row['label']]['x'].append(float(row['paths']))
    series[row['label']]['y'].append(float(row['error']))
    series[row['label']]['slope'] = float(row['slope']) if row['slope'] not in ('', 'nan') else None

def esc(text):
    return str(text).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')

labels = ['MC', 'QMC asset-first', 'QMC time-first', 'QMC PCA+bridge', 'QMC adaptive']
all_x = [x for s in labels for x in series[s]['x']]
all_y = [y for s in labels for y in series[s]['y']]
lx = [__import__('math').log10(v) for v in all_x]
ly = [__import__('math').log10(v) for v in all_y]
xmin, xmax = min(lx), max(lx)
ymin, ymax = min(ly), max(ly)
if xmax == xmin: xmax = xmin + 1
if ymax == ymin: ymax = ymin + 1

svg_w, svg_h = 960, 560
left, right, top, bottom = 86, 24, 60, 64
pw, ph = svg_w - left - right, svg_h - top - bottom
px, py = left, top

parts = [
    f'<svg xmlns="http://www.w3.org/2000/svg" width="{svg_w}" height="{svg_h}" viewBox="0 0 {svg_w} {svg_h}">',
    '<style>',
    'text { font-family: Arial, sans-serif; fill: #111; }',
    '.title { font-size: 18px; font-weight: bold; }',
    '.subtitle { font-size: 12px; fill: #555; }',
    '.axis { font-size: 11px; fill: #444; }',
    '.legend { font-size: 11px; }',
    '.grid { stroke: #cccccc; stroke-width: 1; opacity: 0.45; }',
    '.frame { stroke: #444; stroke-width: 1; fill: none; }',
    '</style>',
    '<text x="480" y="24" text-anchor="middle" class="title">Multi-Asset QMC Dimension Ordering Convergence</text>',
    '<text x="480" y="42" text-anchor="middle" class="subtitle">2-asset spread option under correlated GBM, log-log error vs paths</text>',
]
for i in range(5):
    xv = xmin + (xmax - xmin) * i / 4.0
    xx = px + (xv - xmin) / (xmax - xmin) * pw
    parts.append(f'<line x1="{xx:.2f}" y1="{py}" x2="{xx:.2f}" y2="{py + ph}" class="grid"/>')
    parts.append(f'<text x="{xx:.2f}" y="{py + ph + 18}" text-anchor="middle" class="axis">10^{xv:.2g}</text>')
for i in range(5):
    yv = ymin + (ymax - ymin) * i / 4.0
    yy = py + ph - (yv - ymin) / (ymax - ymin) * ph
    parts.append(f'<line x1="{px}" y1="{yy:.2f}" x2="{px + pw}" y2="{yy:.2f}" class="grid"/>')
    parts.append(f'<text x="{px - 8}" y="{yy + 4:.2f}" text-anchor="end" class="axis">10^{yv:.2g}</text>')
parts.append(f'<rect x="{px}" y="{py}" width="{pw}" height="{ph}" class="frame"/>')

for label in labels:
    pts = []
    for x, y in zip(series[label]['x'], series[label]['y']):
        xx = px + (__import__('math').log10(x) - xmin) / (xmax - xmin) * pw
        yy = py + ph - (__import__('math').log10(y) - ymin) / (ymax - ymin) * ph
        pts.append((xx, yy))
    point_str = ' '.join(f'{xx:.2f},{yy:.2f}' for xx, yy in pts)
    parts.append(f'<polyline points="{point_str}" fill="none" stroke="{colors[label]}" stroke-width="2.2"/>')
    for xx, yy in pts:
        parts.append(f'<circle cx="{xx:.2f}" cy="{yy:.2f}" r="3" fill="{colors[label]}"/>')

guide_specs = [(-0.5, 'MC guide N^-1/2', '#888888'), (-1.0, 'QMC guide N^-1', '#666666')]
anchor_x = max(all_x)
anchor_y = min(all_y)
for slope, guide_label, guide_color in guide_specs:
    guide_pts = []
    c = anchor_y / (anchor_x ** slope)
    for x in sorted(set(all_x)):
        y = c * (x ** slope)
        xx = px + (__import__('math').log10(x) - xmin) / (xmax - xmin) * pw
        yy = py + ph - (__import__('math').log10(y) - ymin) / (ymax - ymin) * ph
        guide_pts.append((xx, yy))
    point_str = ' '.join(f'{xx:.2f},{yy:.2f}' for xx, yy in guide_pts)
    parts.append(f'<polyline points="{point_str}" fill="none" stroke="{guide_color}" stroke-width="1.4" stroke-dasharray="6,5"/>')

parts.append(f'<text x="{svg_w/2}" y="{svg_h - 16}" text-anchor="middle" class="axis">Paths N</text>')
parts.append(f'<text x="20" y="{svg_h/2}" text-anchor="middle" transform="rotate(-90 20 {svg_h/2})" class="axis">Absolute / RMS error</text>')
lx, ly0 = 110, 90
for i, label in enumerate(labels):
    yy = ly0 + i * 16
    slope = series[label]['slope']
    suffix = f' (slope={slope:.2f})' if slope is not None else ''
    parts.append(f'<line x1="{lx}" y1="{yy}" x2="{lx + 18}" y2="{yy}" stroke="{colors[label]}" stroke-width="2.2"/>')
    parts.append(f'<circle cx="{lx + 9}" cy="{yy}" r="3" fill="{colors[label]}"/>')
    parts.append(f'<text x="{lx + 24}" y="{yy + 4}" class="legend">{esc(label + suffix)}</text>')
parts.append('</svg>')
out_svg.write_text('\n'.join(parts), encoding='utf-8')
print(f'Wrote {out_svg}')
)PY";
}

} // namespace

int main() {
    ensureOutputDir();
    writePlotScript();

    const std::vector<double> mu = {0.05, 0.05};
    const std::vector<double> sigma = {0.20, 0.30};
    const double rho = 0.60;
    const std::vector<std::vector<double>> corr = {
        {1.0, rho},
        {rho, 1.0}
    };

    const std::vector<double> s0 = {100.0, 100.0};
    const double r = 0.05;
    const double T = 1.0;
    const double dt = 0.01;
    const double exact = MultiAssetAnalytical::margrabePrice(s0[0], s0[1], sigma[0], sigma[1], rho, T);

    CorrelatedGBM process(mu, sigma, corr);
    SpreadCallPayoff spread(0.0);

    std::vector<SeriesPoint> rows;
    rows.reserve(PATH_COUNTS.size() * 5);

    for (int paths : PATH_COUNTS) {
        rows.push_back(runMcPoint(process, s0, spread, r, T, dt, paths, exact));
        rows.push_back(runQmcPoint("QMC asset-first", process, s0, spread, r, T, dt, paths, exact,
            MultiAssetQMCSimulator::Cholesky, MultiAssetQMCSimulator::AssetFirst));
        rows.push_back(runQmcPoint("QMC time-first", process, s0, spread, r, T, dt, paths, exact,
            MultiAssetQMCSimulator::Cholesky, MultiAssetQMCSimulator::TimeFirst));
        rows.push_back(runQmcPoint("QMC PCA+bridge", process, s0, spread, r, T, dt, paths, exact,
            MultiAssetQMCSimulator::PcaBrownianBridgeHybrid, MultiAssetQMCSimulator::AssetFirst));
        rows.push_back(runQmcPoint("QMC adaptive", process, s0, spread, r, T, dt, paths, exact,
            MultiAssetQMCSimulator::Cholesky, MultiAssetQMCSimulator::AdaptiveVariance));
    }

    for (const std::string label : {
            std::string("MC"),
            std::string("QMC asset-first"),
            std::string("QMC time-first"),
            std::string("QMC PCA+bridge"),
            std::string("QMC adaptive") }) {
        std::vector<int> xs;
        std::vector<double> ys;
        for (const auto& row : rows) {
            if (row.label == label) {
                xs.push_back(row.paths);
                ys.push_back(row.error);
            }
        }
        const double slope = slopeFromLogLog(xs, ys);
        for (auto& row : rows)
            if (row.label == label)
                row.slope = slope;
    }

    std::ofstream csv(CSV_FILE.c_str());
    csv << "label,paths,error,error_stddev,slope\n";
    for (const auto& row : rows) {
        csv << row.label << ","
            << row.paths << ","
            << row.error << ","
            << row.errorStdDev << ","
            << (std::isfinite(row.slope) ? std::to_string(row.slope) : std::string("nan"))
            << "\n";
    }
    csv.close();

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "=== Multi-Asset Ordering Convergence Study ===\n";
    std::cout << "Reference (Margrabe): " << exact << "\n";
    for (const std::string label : {
            std::string("MC"),
            std::string("QMC asset-first"),
            std::string("QMC time-first"),
            std::string("QMC PCA+bridge"),
            std::string("QMC adaptive") }) {
        std::vector<int> xs;
        std::vector<double> ys;
        for (const auto& row : rows) {
            if (row.label == label) {
                xs.push_back(row.paths);
                ys.push_back(row.error);
            }
        }
        std::cout << "  " << label << " slope=" << slopeFromLogLog(xs, ys) << "\n";
        for (const auto& row : rows) {
            if (row.label == label) {
                std::cout << "    N=" << row.paths
                          << " error=" << row.error
                          << " stderr=" << row.errorStdDev
                          << "\n";
            }
        }
    }
    std::cout << "\nWrote: " << CSV_FILE
              << "\nRun: python " << PLOT_SCRIPT << "\n";
    return 0;
}
