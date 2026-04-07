#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

#include "../src/greeks/GreekEstimatorStudy.hpp"
#include "../src/pricing/BlackScholes.hpp"

namespace {

const std::string OUTPUT_DIR = "tests/output";
const std::string CSV_FILE = OUTPUT_DIR + "/qmc_greek_estimators.csv";
const std::string PLOT_SCRIPT = OUTPUT_DIR + "/plot_qmc_greek_estimators.py";

struct SamplerConfig {
    GreekEstimatorStudy::Sampler sampler;
    std::string label;
    bool useOwenScrambling;
};

double slopeFromLogLog(const std::vector<int>& x, const std::vector<double>& y) {
    std::vector<double> xs;
    std::vector<double> ys;
    for (size_t i = 0; i < x.size(); ++i) {
        if (y[i] > 0.0 && std::isfinite(y[i])) {
            xs.push_back(std::log(static_cast<double>(x[i])));
            ys.push_back(std::log(y[i]));
        }
    }
    if (xs.size() < 2)
        return std::numeric_limits<double>::quiet_NaN();

    const double meanX = std::accumulate(xs.begin(), xs.end(), 0.0) / xs.size();
    const double meanY = std::accumulate(ys.begin(), ys.end(), 0.0) / ys.size();

    double num = 0.0;
    double den = 0.0;
    for (size_t i = 0; i < xs.size(); ++i) {
        num += (xs[i] - meanX) * (ys[i] - meanY);
        den += (xs[i] - meanX) * (xs[i] - meanX);
    }
    return num / den;
}

std::string fmt(double x) {
    if (!std::isfinite(x)) return "n/a";
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(3) << x;
    return oss.str();
}

void writePlotScript() {
    std::ofstream py(PLOT_SCRIPT);
    py << "import csv\n";
    py << "import math\n";
    py << "from pathlib import Path\n\n";
    py << "csv_path = Path(r'" << CSV_FILE << "')\n";
    py << R"PY(
out_svg = csv_path.with_suffix('.svg')
rows = []
with csv_path.open(newline='') as f:
    reader = csv.DictReader(f)
    for row in reader:
        if not row or row.get('N') in (None, ''):
            continue
        row['N'] = int(row['N'])
        for key in ['delta_rmse', 'vega_rmse', 'delta_stability', 'vega_stability']:
            row[key] = float(row[key]) if row.get(key) not in (None, '', 'nan') else float('nan')
        for key in ['delta_slope', 'vega_slope']:
            row[key] = float(row[key]) if row.get(key) not in (None, '', 'nan') else float('nan')
        rows.append(row)

metrics = [
    ('delta_rmse', 'Delta RMSE vs Black-Scholes'),
    ('vega_rmse', 'Vega RMSE vs Black-Scholes'),
    ('delta_stability', 'Delta Stability (cross-run stddev)'),
    ('vega_stability', 'Vega Stability (cross-run stddev)'),
]
estimators = ['pathwise', 'likelihood_ratio', 'aad']
sampler_colors = {
    'MC': '#222222',
    'QMC': '#1f77b4',
    'RQMC (Owen)': '#d62728',
}
dash_map = {
    'pathwise': '',
    'likelihood_ratio': '6,4',
    'aad': '2,4',
}
label_map = {
    'pathwise': 'Pathwise',
    'likelihood_ratio': 'Likelihood Ratio',
    'aad': 'AAD',
}

panel_w = 540
panel_h = 280
outer = 22
left = 72
right = 18
top = 38
bottom = 56
plot_w = panel_w - left - right
plot_h = panel_h - top - bottom
svg_w = outer * 3 + panel_w * 2
svg_h = outer * 3 + panel_h * 2

def esc(text):
    return (str(text).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;'))

parts = [
    f'<svg xmlns="http://www.w3.org/2000/svg" width="{svg_w}" height="{svg_h}" viewBox="0 0 {svg_w} {svg_h}">',
    '<style>',
    'text { font-family: Arial, sans-serif; fill: #111; }',
    '.title { font-size: 20px; font-weight: bold; }',
    '.axis { font-size: 11px; fill: #444; }',
    '.panel-title { font-size: 13px; font-weight: bold; }',
    '.legend { font-size: 11px; }',
    '.grid { stroke: #cccccc; stroke-width: 1; opacity: 0.45; }',
    '.frame { stroke: #444; stroke-width: 1; fill: none; }',
    '.series { stroke-width: 2.0; fill: none; }',
    '</style>',
]
parts.append(f'<text x="{svg_w / 2}" y="26" text-anchor="middle" class="title">Pathwise vs Likelihood Ratio vs AAD under MC / QMC / RQMC</text>')

all_n = sorted({row['N'] for row in rows})
for idx, (metric, title) in enumerate(metrics):
    col = idx % 2
    row_idx = idx // 2
    x0 = outer + col * (panel_w + outer)
    y0 = outer + row_idx * (panel_h + outer)
    plot_x0 = x0 + left
    plot_y0 = y0 + top
    panel_rows = [r for r in rows if math.isfinite(r[metric])]
    if not panel_rows:
        continue
    vals = [r[metric] for r in panel_rows if r[metric] > 0]
    y_min = 10 ** math.floor(math.log10(min(vals)))
    y_max = 10 ** math.ceil(math.log10(max(vals)))
    if y_min == y_max:
        y_min /= 10.0
        y_max *= 10.0

    def x_map(x):
        return plot_x0 + (math.log10(x) - math.log10(min(all_n))) / (math.log10(max(all_n)) - math.log10(min(all_n))) * plot_w
    def y_map(y):
        return plot_y0 + plot_h - (math.log10(y) - math.log10(y_min)) / (math.log10(y_max) - math.log10(y_min)) * plot_h

    parts.append(f'<rect x="{x0}" y="{y0}" width="{panel_w}" height="{panel_h}" fill="#ffffff" stroke="#dddddd"/>')
    parts.append(f'<text x="{x0 + panel_w / 2}" y="{y0 + 20}" text-anchor="middle" class="panel-title">{esc(title)}</text>')

    for decade in range(int(math.log10(y_min)), int(math.log10(y_max)) + 1):
        yy = y_map(10 ** decade)
        parts.append(f'<line x1="{plot_x0}" y1="{yy:.2f}" x2="{plot_x0 + plot_w}" y2="{yy:.2f}" class="grid"/>')
        parts.append(f'<text x="{plot_x0 - 8}" y="{yy + 4:.2f}" text-anchor="end" class="axis">1e{decade}</text>')

    for n in all_n:
        xx = x_map(n)
        parts.append(f'<line x1="{xx:.2f}" y1="{plot_y0}" x2="{xx:.2f}" y2="{plot_y0 + plot_h}" class="grid"/>')
        parts.append(f'<text x="{xx:.2f}" y="{plot_y0 + plot_h + 18}" text-anchor="middle" class="axis">{n}</text>')

    parts.append(f'<rect x="{plot_x0}" y="{plot_y0}" width="{plot_w}" height="{plot_h}" class="frame"/>')
    parts.append(f'<text x="{x0 + panel_w / 2}" y="{y0 + panel_h - 14}" text-anchor="middle" class="axis">N (paths)</text>')

    legend_y = plot_y0 + 16
    legend_x = plot_x0 + 8
    legend_i = 0
    for sampler in ['MC', 'QMC', 'RQMC (Owen)']:
        for estimator in estimators:
            series = [r for r in rows if r['sampler'] == sampler and r['estimator'] == estimator and math.isfinite(r[metric]) and r[metric] > 0]
            if not series:
                continue
            series = sorted(series, key=lambda r: r['N'])
            pts = []
            for item in series:
                pts.append((x_map(item['N']), y_map(item[metric])))
            point_str = ' '.join(f'{xx:.2f},{yy:.2f}' for xx, yy in pts)
            dash = dash_map[estimator]
            dash_attr = f' stroke-dasharray="{dash}"' if dash else ''
            color = sampler_colors[sampler]
            parts.append(f'<polyline points="{point_str}" class="series" stroke="{color}"{dash_attr}/>')
            for xx, yy in pts:
                parts.append(f'<circle cx="{xx:.2f}" cy="{yy:.2f}" r="3" fill="{color}"/>')

            ly = legend_y + legend_i * 15
            parts.append(f'<line x1="{legend_x}" y1="{ly}" x2="{legend_x + 18}" y2="{ly}" stroke="{color}" stroke-width="2"{dash_attr}/>')
            slope_vals = [r[metric.replace("rmse", "slope")] for r in rows if r['sampler'] == sampler and r['estimator'] == estimator and 'rmse' in metric and math.isfinite(r.get(metric.replace("rmse", "slope"), float('nan')))]
            slope_txt = f' slope={slope_vals[0]:.2f}' if slope_vals else ''
            parts.append(f'<text x="{legend_x + 24}" y="{ly + 4}" class="legend">{esc(sampler + " / " + label_map[estimator] + slope_txt)}</text>')
            legend_i += 1

parts.append('</svg>')
out_svg.write_text('\n'.join(parts), encoding='utf-8')
print(f'Wrote {out_svg}')
)PY";
}

} // namespace

int main() {
#ifdef _WIN32
    _mkdir(OUTPUT_DIR.c_str());
#else
    mkdir(OUTPUT_DIR.c_str(), 0755);
#endif

    writePlotScript();

    const double S0 = 100.0;
    const double K = 100.0;
    const double r = 0.05;
    const double sigma = 0.20;
    const double T = 1.0;
    const double dt = 0.01;
    const int outerRuns = 8;
    const std::vector<int> pathCounts = {512, 1024, 2048, 4096};
    const Greeks bs = BlackScholes::greeks(S0, K, r, sigma, T);

    const std::vector<SamplerConfig> samplers = {
        {GreekEstimatorStudy::MC, "MC", false},
        {GreekEstimatorStudy::QMC, "QMC", false},
        {GreekEstimatorStudy::RQMC, "RQMC (Owen)", true},
    };

    std::ofstream csv(CSV_FILE);
    csv << "sampler,estimator,N,delta_mean,delta_rmse,delta_stability,delta_inner_stderr,vega_mean,vega_rmse,vega_stability,vega_inner_stderr,delta_slope,vega_slope\n";

    std::cout << "===============================================================\n";
    std::cout << "Greek Estimator Study: Pathwise vs Likelihood Ratio vs AAD\n";
    std::cout << "Sampling schemes: MC, QMC, RQMC\n";
    std::cout << "Targets: Black-Scholes delta and vega\n";
    std::cout << "CSV output: " << CSV_FILE << "\n";
    std::cout << "Plot script: " << PLOT_SCRIPT << "\n";
    std::cout << "===============================================================\n\n";
    std::cout << "Black-Scholes delta = " << std::fixed << std::setprecision(6) << bs.delta
              << ", vega = " << bs.vega << "\n\n";

    for (const auto& sampler : samplers) {
        std::cout << sampler.label << "\n";
        std::cout << std::left
                  << "  " << std::setw(18) << "Estimator"
                  << std::setw(8) << "N"
                  << std::setw(14) << "Delta RMSE"
                  << std::setw(14) << "Delta Stab"
                  << std::setw(14) << "Vega RMSE"
                  << std::setw(14) << "Vega Stab"
                  << std::setw(12) << "dSlope"
                  << std::setw(12) << "vSlope"
                  << "\n";
        std::cout << "  " << std::string(92, '-') << "\n";

        struct NamedSummary {
            std::string name;
            std::vector<double> deltaRmse;
            std::vector<double> vegaRmse;
            std::vector<GreekEstimatorStudy::EstimatorSummary> perN;
        };
        std::vector<NamedSummary> summaries = {
            {"Pathwise", {}, {}, {}},
            {"LikelihoodRatio", {}, {}, {}},
            {"AAD", {}, {}, {}},
        };

        for (int paths : pathCounts) {
            const auto summary = GreekEstimatorStudy::summarize(
                sampler.sampler, K, S0, r, sigma, T, dt, paths, outerRuns,
                "docs/new-joe-kuo-6.21201.txt", true, sampler.useOwenScrambling, 8, 42);

            summaries[0].perN.push_back(summary.pathwise);
            summaries[1].perN.push_back(summary.likelihoodRatio);
            summaries[2].perN.push_back(summary.aad);
            summaries[0].deltaRmse.push_back(summary.pathwise.deltaRmse);
            summaries[1].deltaRmse.push_back(summary.likelihoodRatio.deltaRmse);
            summaries[2].deltaRmse.push_back(summary.aad.deltaRmse);
            summaries[0].vegaRmse.push_back(summary.pathwise.vegaRmse);
            summaries[1].vegaRmse.push_back(summary.likelihoodRatio.vegaRmse);
            summaries[2].vegaRmse.push_back(summary.aad.vegaRmse);
        }

        for (const auto& entry : summaries) {
            const double deltaSlope = slopeFromLogLog(pathCounts, entry.deltaRmse);
            const double vegaSlope = slopeFromLogLog(pathCounts, entry.vegaRmse);
            for (size_t i = 0; i < pathCounts.size(); ++i) {
                const auto& stats = entry.perN[i];
                csv << sampler.label << ","
                    << entry.name << ","
                    << pathCounts[i] << ","
                    << stats.deltaMean << ","
                    << stats.deltaRmse << ","
                    << (std::isfinite(stats.deltaStability) ? std::to_string(stats.deltaStability) : std::string("nan")) << ","
                    << stats.deltaInnerStdErr << ","
                    << stats.vegaMean << ","
                    << stats.vegaRmse << ","
                    << (std::isfinite(stats.vegaStability) ? std::to_string(stats.vegaStability) : std::string("nan")) << ","
                    << stats.vegaInnerStdErr << ","
                    << (std::isfinite(deltaSlope) ? std::to_string(deltaSlope) : std::string("nan")) << ","
                    << (std::isfinite(vegaSlope) ? std::to_string(vegaSlope) : std::string("nan")) << "\n";

                std::cout << "  " << std::setw(18) << entry.name
                          << std::setw(8) << pathCounts[i]
                          << std::setw(14) << fmt(stats.deltaRmse)
                          << std::setw(14) << fmt(stats.deltaStability)
                          << std::setw(14) << fmt(stats.vegaRmse)
                          << std::setw(14) << fmt(stats.vegaStability)
                          << std::setw(12) << (i == 0 ? fmt(deltaSlope) : "")
                          << std::setw(12) << (i == 0 ? fmt(vegaSlope) : "")
                          << "\n";
            }
        }
        std::cout << "\n";
    }

    csv.close();

    const auto mcCheck = GreekEstimatorStudy::summarize(
        GreekEstimatorStudy::MC, K, S0, r, sigma, T, dt, 4096, outerRuns);
    assert(mcCheck.pathwise.deltaRmse < 0.05);
    assert(mcCheck.aad.deltaRmse < 0.05);
    assert(mcCheck.likelihoodRatio.deltaRmse < 0.12);
    assert(mcCheck.pathwise.vegaRmse < 2.0);
    assert(mcCheck.aad.vegaRmse < 2.0);
    assert(mcCheck.likelihoodRatio.vegaRmse < 6.0);

    std::cout << "Finished. Run the plot script with:\n";
    std::cout << "  python " << PLOT_SCRIPT << "\n";
    return 0;
}
