#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <string>
#include <vector>

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

#include "../src/roughvol/RoughVolatilityResearch.hpp"

namespace {

const std::string OUTPUT_DIR = "tests/output";
const std::string PRICING_CSV = OUTPUT_DIR + "/rough_stress_pricing.csv";
const std::string SUMMARY_CSV = OUTPUT_DIR + "/rough_stress_summary.csv";
const std::string PLOT_SCRIPT = OUTPUT_DIR + "/plot_rough_stress_robustness.py";

struct Scenario {
    std::string label;
    std::string family;
    RoughVolatilityResearch::RoughBergomiParameters params;
};

struct PricingPoint {
    std::string scenario;
    std::string family;
    std::string sampler;
    int paths;
    double rmse;
    double rmseStdDev;
};

struct SummaryPoint {
    std::string scenario;
    std::string family;
    double atmIv;
    double leftSkew;
    double empiricalRho;
    double rhoAbsError;
    double mcSlope;
    double qmcSlope;
    double rqmcSlope;
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

double logLogSlope(const std::vector<double>& x, const std::vector<double>& y) {
    std::vector<double> xs;
    std::vector<double> ys;
    for (std::size_t i = 0; i < x.size(); ++i) {
        if (x[i] > 0.0 && y[i] > 0.0 && std::isfinite(y[i])) {
            xs.push_back(std::log(x[i]));
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

void writePlotScript() {
    std::ofstream py(PLOT_SCRIPT.c_str());
    py << "import csv\n";
    py << "from pathlib import Path\n\n";
    py << "base = Path(r'" << OUTPUT_DIR << "')\n";
    py << R"PY(
pricing_rows = list(csv.DictReader((base / 'rough_stress_pricing.csv').open(newline='')))
summary_rows = list(csv.DictReader((base / 'rough_stress_summary.csv').open(newline='')))

scenario_order = [row['scenario'] for row in summary_rows]
family_by_scenario = {row['scenario']: row['family'] for row in summary_rows}
sampler_colors = {'MC': '#222222', 'QMC': '#1f77b4', 'RQMC': '#d62728'}
family_colors = {'baseline': '#9ecae1', 'correlation': '#fdd0a2', 'vol-of-vol': '#fdae6b', 'roughness': '#c7e9c0'}

def esc(text):
    return str(text).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')

def grouped_bars(parts, panel_x, panel_y, panel_w, panel_h, title, values_by_key, y_label, guide_lines=None):
    left, right, top, bottom = 64, 20, 40, 78
    plot_w = panel_w - left - right
    plot_h = panel_h - top - bottom
    x0 = panel_x + left
    y0 = panel_y + top
    values = [v for v in values_by_key.values() if v is not None]
    if not values:
        return
    if guide_lines:
        values.extend(guide_lines)
    ymin = min(values)
    ymax = max(values)
    if ymax == ymin:
        ymax = ymin + 1.0
    pad = 0.08 * (ymax - ymin)
    ymin -= pad
    ymax += pad
    parts.append(f'<rect x="{panel_x}" y="{panel_y}" width="{panel_w}" height="{panel_h}" fill="#fff" stroke="#dddddd"/>')
    parts.append(f'<text x="{panel_x + panel_w/2}" y="{panel_y + 20}" text-anchor="middle" class="legend">{esc(title)}</text>')
    for i in range(5):
        yv = ymin + (ymax - ymin) * i / 4.0
        yy = y0 + plot_h - (yv - ymin) / (ymax - ymin) * plot_h
        parts.append(f'<line x1="{x0}" y1="{yy:.2f}" x2="{x0 + plot_w}" y2="{yy:.2f}" class="grid"/>')
        parts.append(f'<text x="{x0 - 8}" y="{yy + 4:.2f}" text-anchor="end" class="axis">{yv:.3g}</text>')
    if guide_lines:
        for gv in guide_lines:
            yy = y0 + plot_h - (gv - ymin) / (ymax - ymin) * plot_h
            parts.append(f'<line x1="{x0}" y1="{yy:.2f}" x2="{x0 + plot_w}" y2="{yy:.2f}" stroke="#888" stroke-dasharray="5,4" stroke-width="1.2"/>')
    parts.append(f'<rect x="{x0}" y="{y0}" width="{plot_w}" height="{plot_h}" class="frame"/>')
    group_w = plot_w / max(1, len(scenario_order))
    bar_w = group_w * 0.18
    offsets = {'MC': -bar_w * 1.2, 'QMC': 0.0, 'RQMC': bar_w * 1.2}
    for i, scenario in enumerate(scenario_order):
        gx = x0 + (i + 0.5) * group_w
        bg = family_colors.get(family_by_scenario.get(scenario, ''), '#f4f4f4')
        parts.append(f'<rect x="{x0 + i * group_w + 2:.2f}" y="{y0}" width="{group_w - 4:.2f}" height="{plot_h}" fill="{bg}" opacity="0.25"/>')
        for sampler in ['MC', 'QMC', 'RQMC']:
            value = values_by_key.get((scenario, sampler))
            if value is None:
                continue
            yy = y0 + plot_h - (value - ymin) / (ymax - ymin) * plot_h
            base = y0 + plot_h - (0.0 - ymin) / (ymax - ymin) * plot_h if ymin <= 0.0 <= ymax else (y0 + plot_h)
            top_y = min(yy, base)
            height = abs(base - yy)
            parts.append(f'<rect x="{gx + offsets[sampler] - bar_w/2:.2f}" y="{top_y:.2f}" width="{bar_w:.2f}" height="{height:.2f}" fill="{sampler_colors[sampler]}"/>')
        parts.append(f'<text x="{gx:.2f}" y="{y0 + plot_h + 18}" text-anchor="middle" class="axis">{esc(scenario)}</text>')
    parts.append(f'<text x="{panel_x + panel_w/2}" y="{panel_y + panel_h - 14}" text-anchor="middle" class="axis">Stress scenario</text>')
    parts.append(f'<text x="{panel_x + 18}" y="{panel_y + panel_h/2}" text-anchor="middle" transform="rotate(-90 {panel_x + 18} {panel_y + panel_h/2})" class="axis">{esc(y_label)}</text>')

pricing_slopes = {}
final_rmse = {}
for row in pricing_rows:
    pricing_slopes[(row['scenario'], row['sampler'])] = float(row['slope']) if row['slope'] not in ('', 'nan') else None
    key = (row['scenario'], row['sampler'])
    current = final_rmse.get(key)
    if current is None or int(row['paths']) > current[0]:
        final_rmse[key] = (int(row['paths']), float(row['rmse']))
final_rmse = {k: v[1] for k, v in final_rmse.items()}

improvement = {}
severity = {}
for row in summary_rows:
    scenario = row['scenario']
    improvement[(scenario, 'MC')] = 0.0
    improvement[(scenario, 'QMC')] = float(row['qmc_minus_mc_slope'])
    improvement[(scenario, 'RQMC')] = float(row['rqmc_minus_mc_slope'])
    severity[(scenario, 'MC')] = float(row['atm_iv'])
    severity[(scenario, 'QMC')] = float(row['left_skew'])
    severity[(scenario, 'RQMC')] = float(row['rho_abs_error'])

parts = [
    '<svg xmlns="http://www.w3.org/2000/svg" width="1260" height="940" viewBox="0 0 1260 940">',
    '<style>',
    'text { font-family: Arial, sans-serif; fill: #111; }',
    '.title { font-size: 19px; font-weight: bold; }',
    '.subtitle { font-size: 12px; fill: #555; }',
    '.axis { font-size: 11px; fill: #444; }',
    '.legend { font-size: 11px; }',
    '.grid { stroke: #d0d0d0; stroke-width: 1; opacity: 0.5; }',
    '.frame { stroke: #444; stroke-width: 1; fill: none; }',
    '</style>',
    '<text x="630" y="24" text-anchor="middle" class="title">Rough Bergomi Stress Robustness</text>',
    '<text x="630" y="42" text-anchor="middle" class="subtitle">Extreme correlation, high vol-of-vol, and very small H show where QMC drifts away from ideal behavior</text>',
]

grouped_bars(parts, 20, 60, 1220, 270, 'Pricing convergence slope by scenario', pricing_slopes, 'Log-log slope of pricing RMSE vs paths', [-1.0, -0.5])
grouped_bars(parts, 20, 350, 1220, 250, 'Final-N pricing RMSE by scenario', final_rmse, 'RMSE at largest path count')
grouped_bars(parts, 20, 620, 600, 260, 'QMC / RQMC slope advantage over MC', improvement, 'Slope minus MC slope')
grouped_bars(parts, 640, 620, 600, 260, 'Regime severity proxies', severity, 'ATM IV / left skew / |rho error|')

lx, ly = 90, 78
for idx, sampler in enumerate(['MC', 'QMC', 'RQMC']):
    yy = ly + idx * 16
    parts.append(f'<line x1="{lx}" y1="{yy}" x2="{lx + 18}" y2="{yy}" stroke="{sampler_colors[sampler]}" stroke-width="2.2"/>')
    parts.append(f'<text x="{lx + 24}" y="{yy + 4}" class="legend">{sampler}</text>')

fx, fy = 980, 78
families = [('baseline', 'Baseline'), ('correlation', 'Extreme corr'), ('vol-of-vol', 'High vol-of-vol'), ('roughness', 'Very small H')]
for idx, (key, label) in enumerate(families):
    yy = fy + idx * 16
    parts.append(f'<rect x="{fx}" y="{yy - 9}" width="14" height="10" fill="{family_colors[key]}" opacity="0.55"/>')
    parts.append(f'<text x="{fx + 22}" y="{yy}" class="legend">{esc(label)}</text>')

parts.append('</svg>')
(base / 'rough_stress_robustness.svg').write_text('\n'.join(parts), encoding='utf-8')
print('Wrote', base / 'rough_stress_robustness.svg')
)PY";
}

PricingPoint evaluatePricingPoint(
    const Scenario& scenario,
    const char* samplerLabel,
    RoughBergomiModel::SamplingMethod method,
    bool useOwen,
    int paths,
    int outerRuns,
    double benchmarkPrice,
    double strike,
    double maturity,
    double rate,
    int steps)
{
    const int runs = (method == RoughBergomiModel::QMC) ? 1 : outerRuns;
    std::vector<double> errors;
    errors.reserve(runs);

    for (int run = 0; run < runs; ++run) {
        const auto result = RoughBergomiModel::priceEuropeanCall(
            scenario.params,
            100.0,
            rate,
            maturity,
            steps,
            paths,
            strike,
            method,
            "docs/new-joe-kuo-6.21201.txt",
            useOwen,
            8,
            static_cast<std::uint64_t>(1000 + paths * 31 + run * 101));
        errors.push_back(std::abs(result.price - benchmarkPrice));
    }

    PricingPoint point;
    point.scenario = scenario.label;
    point.family = scenario.family;
    point.sampler = samplerLabel;
    point.paths = paths;
    double sq = 0.0;
    for (double e : errors)
        sq += e * e;
    point.rmse = std::sqrt(sq / std::max<std::size_t>(1, errors.size()));
    point.rmseStdDev = stdDev(errors);
    return point;
}

} // namespace

int main(int argc, char** argv) {
    std::cout << std::unitbuf;
    ensureOutputDir();
    writePlotScript();

    const bool fullMode = (argc > 1 && std::string(argv[1]) == "full");
    const double maturity = 1.0;
    const double rate = 0.01;
    const double strike = 100.0;
    const int steps = fullMode ? 12 : 10;
    const int benchmarkPaths = fullMode ? 2048 : 512;
    const std::vector<int> pricingPathCounts = fullMode
        ? std::vector<int>{64, 128, 256, 512}
        : std::vector<int>{64, 128, 256};
    const int pricingOuterRuns = fullMode ? 4 : 2;

    const std::vector<Scenario> scenarios = {
        {"Baseline", "baseline", {0.04, 1.5, -0.70, 0.12}},
        {"Rho=-0.95", "correlation", {0.04, 1.5, -0.95, 0.12}},
        {"Rho=+0.95", "correlation", {0.04, 1.5,  0.95, 0.12}},
        {"Eta=3.0", "vol-of-vol", {0.04, 3.0, -0.70, 0.12}},
        {"H=0.03", "roughness", {0.04, 1.5, -0.70, 0.03}},
    };

    RoughVolatilityResearch::RoughSurfaceSpec surfaceSpec;
    surfaceSpec.spot = 100.0;
    surfaceSpec.rate = rate;
    surfaceSpec.maturity = maturity;
    surfaceSpec.steps = steps;
    surfaceSpec.paths = fullMode ? 192 : 64;
    surfaceSpec.rqmcReplicates = 8;
    surfaceSpec.strikes = {90.0, 100.0, 110.0};
    surfaceSpec.expiries = {0.5, 1.0};

    std::vector<PricingPoint> pricingRows;
    std::vector<SummaryPoint> summaryRows;

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "=== Rough Volatility Stress / Robustness Study ===\n";
    std::cout << "Mode: " << (fullMode ? "full" : "fast") << "\n";

    const struct SamplerSpec {
        const char* label;
        RoughBergomiModel::SamplingMethod method;
        bool useOwen;
    } samplers[] = {
        {"MC", RoughBergomiModel::MC, false},
        {"QMC", RoughBergomiModel::QMC, false},
        {"RQMC", RoughBergomiModel::RQMC, true},
    };

    for (const Scenario& scenario : scenarios) {
        std::cout << "\nScenario: " << scenario.label
                  << "  [xi0=" << scenario.params.xi0
                  << ", eta=" << scenario.params.eta
                  << ", rho=" << scenario.params.rho
                  << ", H=" << scenario.params.hurst << "]\n";
        try {
            const auto benchmark = RoughBergomiModel::priceEuropeanCall(
                scenario.params,
                100.0,
                rate,
                maturity,
                steps,
                benchmarkPaths,
                strike,
                RoughBergomiModel::RQMC,
                "docs/new-joe-kuo-6.21201.txt",
                true,
                8,
                7777ULL);

            std::map<std::string, std::vector<double> > slopeX;
            std::map<std::string, std::vector<double> > slopeY;
            for (const auto& sampler : samplers) {
                for (int paths : pricingPathCounts) {
                    const PricingPoint point = evaluatePricingPoint(
                        scenario,
                        sampler.label,
                        sampler.method,
                        sampler.useOwen,
                        paths,
                        pricingOuterRuns,
                        benchmark.price,
                        strike,
                        maturity,
                        rate,
                        steps);
                    pricingRows.push_back(point);
                    slopeX[sampler.label].push_back(static_cast<double>(paths));
                    slopeY[sampler.label].push_back(std::max(point.rmse, 1e-12));
                }
            }

            const double mcSlope = logLogSlope(slopeX["MC"], slopeY["MC"]);
            const double qmcSlope = logLogSlope(slopeX["QMC"], slopeY["QMC"]);
            const double rqmcSlope = logLogSlope(slopeX["RQMC"], slopeY["RQMC"]);
            std::cout << "  MC slope=" << mcSlope
                      << " QMC slope=" << qmcSlope
                      << " RQMC slope=" << rqmcSlope << "\n";

            const auto corr = RoughBergomiModel::empiricalDriverCorrelation(
                scenario.params, maturity, steps, 200, 4321ULL);

            const VolSurface target = RoughVolatilityResearch::generateRoughBergomiSurface(
                scenario.params, surfaceSpec, 123, RoughBergomiModel::RQMC, true);
            const auto smile = RoughVolatilityResearch::smileSlice(target, 1.0, surfaceSpec.strikes);
            const double leftSkew = RoughVolatilityResearch::smileLeftSkewMetric(smile);
            const double atmIv = target.impliedVol(1.0, 100.0);

            std::cout << "  atmIV=" << atmIv
                      << " leftSkew=" << leftSkew
                      << " empirical rho=" << corr.empiricalRho << "\n";

            summaryRows.push_back({
                scenario.label,
                scenario.family,
                atmIv,
                leftSkew,
                corr.empiricalRho,
                std::abs(corr.empiricalRho - corr.targetRho),
                mcSlope,
                qmcSlope,
                rqmcSlope
            });
        } catch (const std::exception& e) {
            std::cout << "  scenario failed: " << e.what() << "\n";
        }
    }

    std::ofstream pricing(PRICING_CSV.c_str());
    pricing << "scenario,family,sampler,paths,rmse,rmse_stddev,slope\n";
    for (const auto& scenario : scenarios) {
        for (const char* sampler : {"MC", "QMC", "RQMC"}) {
            std::vector<double> xs;
            std::vector<double> ys;
            for (const auto& row : pricingRows) {
                if (row.scenario == scenario.label && row.sampler == sampler) {
                    xs.push_back(static_cast<double>(row.paths));
                    ys.push_back(row.rmse);
                }
            }
            const double slope = logLogSlope(xs, ys);
            for (const auto& row : pricingRows) {
                if (row.scenario == scenario.label && row.sampler == sampler) {
                    pricing << row.scenario << ","
                            << row.family << ","
                            << row.sampler << ","
                            << row.paths << ","
                            << row.rmse << ","
                            << row.rmseStdDev << ","
                            << (std::isfinite(slope) ? std::to_string(slope) : std::string("nan"))
                            << "\n";
                }
            }
        }
    }
    pricing.close();

    std::ofstream summary(SUMMARY_CSV.c_str());
    summary << "scenario,family,atm_iv,left_skew,empirical_rho,rho_abs_error,mc_slope,qmc_slope,rqmc_slope,qmc_minus_mc_slope,rqmc_minus_mc_slope\n";
    for (const auto& row : summaryRows) {
        summary << row.scenario << ","
                << row.family << ","
                << row.atmIv << ","
                << row.leftSkew << ","
                << row.empiricalRho << ","
                << row.rhoAbsError << ","
                << row.mcSlope << ","
                << row.qmcSlope << ","
                << row.rqmcSlope << ","
                << (row.qmcSlope - row.mcSlope) << ","
                << (row.rqmcSlope - row.mcSlope) << "\n";
    }
    summary.close();

    std::cout << "\nWrote:\n  " << PRICING_CSV
              << "\n  " << SUMMARY_CSV
              << "\nRun: python " << PLOT_SCRIPT
              << "\nOptional heavier run: .\\test_rough_stress_robustness.exe full\n";
    return 0;
}
