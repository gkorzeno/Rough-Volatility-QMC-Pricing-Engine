#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "../src/core/Random.hpp"
#include "../src/roughvol/RoughVolatilityResearch.hpp"
#include "../src/surface/SVI.hpp"

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

namespace {

const std::string OUTPUT_DIR = "tests/output";

void ensureOutputDir() {
#ifdef _WIN32
    _mkdir(OUTPUT_DIR.c_str());
#else
    mkdir(OUTPUT_DIR.c_str(), 0755);
#endif
}

void writePlotScript() {
    std::ofstream py((OUTPUT_DIR + "/plot_roughvol_diagnostics.py").c_str());
    py << "import csv\n";
    py << "from pathlib import Path\n\n";
    py << "base = Path(r'" << OUTPUT_DIR << "')\n";
    py << R"PY(

def read_rows(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f))

def esc(text):
    return str(text).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')

def line_points(xs, ys, x0, y0, w, h):
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    if xmax == xmin:
        xmax = xmin + 1.0
    if ymax == ymin:
        ymax = ymin + 1.0e-6
    pts = []
    for x, y in zip(xs, ys):
        xx = x0 + (x - xmin) / (xmax - xmin) * w
        yy = y0 + h - (y - ymin) / (ymax - ymin) * h
        pts.append((xx, yy))
    return pts, xmin, xmax, ymin, ymax

def render_panel(title, xlabel, ylabel, series, out_path, width=700, height=420):
    margin_l = 70
    margin_r = 20
    margin_t = 45
    margin_b = 55
    plot_w = width - margin_l - margin_r
    plot_h = height - margin_t - margin_b
    x0 = margin_l
    y0 = margin_t
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<style>',
        'text { font-family: Arial, sans-serif; fill: #111; }',
        '.title { font-size: 18px; font-weight: bold; }',
        '.axis { font-size: 11px; fill: #444; }',
        '.legend { font-size: 11px; }',
        '.grid { stroke: #cccccc; stroke-width: 1; opacity: 0.45; }',
        '.frame { stroke: #444; stroke-width: 1; fill: none; }',
        '</style>',
        f'<text x="{width/2}" y="24" text-anchor="middle" class="title">{esc(title)}</text>',
    ]
    all_x = []
    all_y = []
    computed = []
    for label, color, xs, ys in series:
        pts, xmin, xmax, ymin, ymax = line_points(xs, ys, x0, y0, plot_w, plot_h)
        computed.append((label, color, pts))
        all_x.extend(xs)
        all_y.extend(ys)
    xmin, xmax = min(all_x), max(all_x)
    ymin, ymax = min(all_y), max(all_y)
    if xmax == xmin:
        xmax = xmin + 1.0
    if ymax == ymin:
        ymax = ymin + 1.0e-6
    for i in range(5):
        xv = xmin + (xmax - xmin) * i / 4.0
        xx = x0 + (xv - xmin) / (xmax - xmin) * plot_w
        parts.append(f'<line x1="{xx:.2f}" y1="{y0}" x2="{xx:.2f}" y2="{y0 + plot_h}" class="grid"/>')
        parts.append(f'<text x="{xx:.2f}" y="{y0 + plot_h + 18}" text-anchor="middle" class="axis">{xv:.3g}</text>')
    for i in range(5):
        yv = ymin + (ymax - ymin) * i / 4.0
        yy = y0 + plot_h - (yv - ymin) / (ymax - ymin) * plot_h
        parts.append(f'<line x1="{x0}" y1="{yy:.2f}" x2="{x0 + plot_w}" y2="{yy:.2f}" class="grid"/>')
        parts.append(f'<text x="{x0 - 8}" y="{yy + 4:.2f}" text-anchor="end" class="axis">{yv:.3g}</text>')
    parts.append(f'<rect x="{x0}" y="{y0}" width="{plot_w}" height="{plot_h}" class="frame"/>')
    for label, color, pts in computed:
        point_str = ' '.join(f'{xx:.2f},{yy:.2f}' for xx, yy in pts)
        parts.append(f'<polyline points="{point_str}" fill="none" stroke="{color}" stroke-width="2.2"/>')
        for xx, yy in pts:
            parts.append(f'<circle cx="{xx:.2f}" cy="{yy:.2f}" r="3" fill="{color}"/>')
    parts.append(f'<text x="{width/2}" y="{height - 12}" text-anchor="middle" class="axis">{esc(xlabel)}</text>')
    parts.append(f'<text x="18" y="{height/2}" text-anchor="middle" transform="rotate(-90 18 {height/2})" class="axis">{esc(ylabel)}</text>')
    lx = x0 + 8
    ly = y0 + 14
    for idx, item in enumerate(computed):
        label, color, _ = item
        yy = ly + idx * 16
        parts.append(f'<line x1="{lx}" y1="{yy}" x2="{lx + 18}" y2="{yy}" stroke="{color}" stroke-width="2.2"/>')
        parts.append(f'<circle cx="{lx + 9}" cy="{yy}" r="3" fill="{color}"/>')
        parts.append(f'<text x="{lx + 24}" y="{yy + 4}" class="legend">{esc(label)}</text>')
    parts.append('</svg>')
    out_path.write_text('\n'.join(parts), encoding='utf-8')

smile = read_rows(base / 'rough_smile.csv')
term = read_rows(base / 'rough_term_structure.csv')
svi = read_rows(base / 'rough_vs_svi.csv')
var = read_rows(base / 'rough_mc_variance.csv')
calib = read_rows(base / 'rough_calibration_noise.csv')

render_panel(
    'Rough Bergomi Smile Diagnostics',
    'Strike',
    'Implied Vol',
    [('Rough Smile', '#1f77b4',
      [float(r['strike']) for r in smile],
      [float(r['implied_vol']) for r in smile])],
    base / 'rough_smile.svg'
)

render_panel(
    'ATM Term Structure',
    'Maturity',
    'Implied Vol',
    [('ATM IV', '#d62728',
      [float(r['expiry']) for r in term],
      [float(r['implied_vol']) for r in term])],
    base / 'rough_term_structure.svg'
)

svi_groups = {}
for row in svi:
    key = row['expiry']
    svi_groups.setdefault(key, {'x': [], 'rough': [], 'svi': []})
    svi_groups[key]['x'].append(float(row['strike']))
    svi_groups[key]['rough'].append(float(row['rough_iv']))
    svi_groups[key]['svi'].append(float(row['svi_iv']))
first_key = sorted(svi_groups.keys())[0]
render_panel(
    f'Rough Surface vs SVI Fit (T={first_key})',
    'Strike',
    'Implied Vol',
    [
        ('Rough', '#1f77b4', svi_groups[first_key]['x'], svi_groups[first_key]['rough']),
        ('SVI', '#2ca02c', svi_groups[first_key]['x'], svi_groups[first_key]['svi']),
    ],
    base / 'rough_vs_svi.svg'
)

render_panel(
    'Monte Carlo Variance Analysis',
    'Paths',
    'Surface RMSE',
    [('RMSE', '#9467bd',
      [float(r['paths']) for r in var],
      [float(r['rmse']) for r in var])],
    base / 'rough_mc_variance.svg'
)

if calib:
    by_sampler = {}
    for row in calib:
        by_sampler.setdefault(row['sampler'], {'x': [], 'prmse': [], 'ormse': []})
        by_sampler[row['sampler']]['x'].append(float(row['paths']))
        by_sampler[row['sampler']]['prmse'].append(float(row['parameter_rmse_mean']))
        by_sampler[row['sampler']]['ormse'].append(float(row['objective_rmse_mean']))

    render_panel(
        'Calibration Under Noise: Parameter RMSE',
        'Paths per evaluation',
        'Parameter RMSE',
        [(label, color, vals['x'], vals['prmse']) for label, color, vals in [
            ('MC', '#222222', by_sampler.get('MC', {'x': [], 'prmse': []})),
            ('QMC', '#1f77b4', by_sampler.get('QMC', {'x': [], 'prmse': []})),
            ('RQMC', '#d62728', by_sampler.get('RQMC', {'x': [], 'prmse': []})),
        ] if vals['x']],
        base / 'rough_calibration_noise.svg'
    )

print('Wrote rough-vol SVG diagnostics to', base)
)PY";
}

double logLogSlope(const std::vector<RoughVolatilityResearch::VariancePoint>& variance) {
    if (variance.size() < 2)
        return 0.0;
    double mx = 0.0, my = 0.0;
    for (std::size_t i = 0; i < variance.size(); ++i) {
        mx += std::log(static_cast<double>(variance[i].paths));
        my += std::log(std::max(variance[i].rmse, 1e-12));
    }
    mx /= variance.size();
    my /= variance.size();
    double num = 0.0, den = 0.0;
    for (std::size_t i = 0; i < variance.size(); ++i) {
        const double x = std::log(static_cast<double>(variance[i].paths));
        const double y = std::log(std::max(variance[i].rmse, 1e-12));
        num += (x - mx) * (y - my);
        den += (x - mx) * (x - mx);
    }
    return den > 0.0 ? num / den : 0.0;
}

} // namespace

int main() {
    std::cout << std::unitbuf;
    ensureOutputDir();
    writePlotScript();

    Random rng(42);

    const double T = 1.0;
    const int steps = 12;
    const double hurst = 0.12;

    auto fBm = RoughVolatilityResearch::sampleFractionalBrownianMotion(T, steps, hurst, rng);
    auto fractionalPath = RoughVolatilityResearch::simulateFractionalSDE(
        1.0,
        T,
        fBm,
        [](double x, double) { return -0.3 * x; },
        [](double, double) { return 0.5; });

    RoughVolatilityResearch::RoughBergomiParameters params;
    params.xi0 = 0.04;
    params.eta = 1.6;
    params.rho = -0.75;
    params.hurst = hurst;

    RoughVolatilityResearch::RoughSurfaceSpec surfaceSpec;
    surfaceSpec.spot = 100.0;
    surfaceSpec.rate = 0.01;
    surfaceSpec.maturity = T;
    surfaceSpec.steps = steps;
    surfaceSpec.paths = 40;
    surfaceSpec.rqmcReplicates = 8;
    surfaceSpec.strikes = {90.0, 100.0, 110.0};
    surfaceSpec.expiries = {0.50, 1.00};

    std::cout << "[progress] generating target rough surface\n";
    VolSurface target = RoughVolatilityResearch::generateRoughBergomiSurface(params, surfaceSpec, 42);

    RoughVolatilityResearch::CalibrationSpec calibSpec;
    calibSpec.spot = surfaceSpec.spot;
    calibSpec.rate = surfaceSpec.rate;
    calibSpec.maturity = surfaceSpec.maturity;
    calibSpec.steps = surfaceSpec.steps;
    calibSpec.pathsPerEvaluation = 16;
    calibSpec.maxIterations = 4;
    calibSpec.tolerance = 1e-4;
    calibSpec.objectiveType = RoughBergomiCalibrator::SviSmoothedSurfaceRMSE;
    calibSpec.initialGuess = {0.05, 1.3, -0.50, 0.18};
    calibSpec.lowerBounds = {0.01, 0.20, -0.99, 0.02};
    calibSpec.upperBounds = {0.20, 3.00, -0.05, 0.49};

    std::cout << "[progress] calibrating rough Bergomi\n";
    auto fit = RoughVolatilityResearch::calibrateRoughBergomi(target, calibSpec);
    std::cout << "[progress] computing smile and term-structure slices\n";
    auto smile = RoughVolatilityResearch::smileSlice(target, 1.0, surfaceSpec.strikes);
    auto term = RoughVolatilityResearch::termStructureSlice(target, surfaceSpec.spot);
    double leftSkew = RoughVolatilityResearch::smileLeftSkewMetric(smile);
    std::cout << "[progress] fitting SVI slice at T=1.0\n";
    const double sviExpiry = 1.0;
    const double F = surfaceSpec.spot * std::exp(surfaceSpec.rate * sviExpiry);
    std::vector<double> marketVols;
    for (std::size_t i = 0; i < surfaceSpec.strikes.size(); ++i)
        marketVols.push_back(target.impliedVol(sviExpiry, surfaceSpec.strikes[i]));
    SVIParams sviParams = SVI::fit(surfaceSpec.strikes, marketVols, F, sviExpiry);
    double sviRmse = 0.0;
    double sviMaxErr = 0.0;
    {
        std::ofstream out((OUTPUT_DIR + "/rough_vs_svi.csv").c_str());
        out << "expiry,strike,rough_iv,svi_iv,error\n";
        for (std::size_t i = 0; i < surfaceSpec.strikes.size(); ++i) {
            const double roughIv = target.impliedVol(sviExpiry, surfaceSpec.strikes[i]);
            double sviIv = roughIv;
            try {
                sviIv = SVI::impliedVol(sviParams, surfaceSpec.strikes[i], F, sviExpiry);
            } catch (...) {
                sviIv = roughIv;
            }
            const double err = sviIv - roughIv;
            sviRmse += err * err;
            sviMaxErr = std::max(sviMaxErr, std::abs(err));
            out << sviExpiry << "," << surfaceSpec.strikes[i] << ","
                << roughIv << "," << sviIv << "," << err << "\n";
        }
        sviRmse = std::sqrt(sviRmse / surfaceSpec.strikes.size());
    }
    std::cout << "[progress] running variance analysis\n";
    auto variance = RoughVolatilityResearch::varianceAnalysis(
        params, surfaceSpec, std::vector<int>{20, 40, 80, 160}, 400, 4);
    std::cout << "[progress] pricing rough Bergomi under MC/QMC/RQMC\n";
    const double pricingStrike = 100.0;
    RoughBergomiModel::PricingResult mcPrice = RoughBergomiModel::priceEuropeanCall(
        params, surfaceSpec.spot, surfaceSpec.rate, 1.0, steps, 512, pricingStrike,
        RoughBergomiModel::MC);
    RoughBergomiModel::PricingResult qmcPrice = RoughBergomiModel::priceEuropeanCall(
        params, surfaceSpec.spot, surfaceSpec.rate, 1.0, steps, 512, pricingStrike,
        RoughBergomiModel::QMC);
    RoughBergomiModel::PricingResult rqmcPrice = RoughBergomiModel::priceEuropeanCall(
        params, surfaceSpec.spot, surfaceSpec.rate, 1.0, steps, 512, pricingStrike,
        RoughBergomiModel::RQMC, "docs/new-joe-kuo-6.21201.txt", true, 8);
    std::cout << "[progress] checking correlated drivers\n";
    RoughBergomiModel::DriverCorrelationResult corr =
        RoughBergomiModel::empiricalDriverCorrelation(params, 1.0, steps, 400, 321ULL);

    std::cout << "[progress] writing diagnostics\n";
    RoughVolDiagnostics::writeSmileCsv(OUTPUT_DIR + "/rough_smile.csv", smile);
    RoughVolDiagnostics::writeTermCsv(OUTPUT_DIR + "/rough_term_structure.csv", term);
    RoughVolDiagnostics::writeVarianceCsv(OUTPUT_DIR + "/rough_mc_variance.csv", variance);

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "=== Rough Volatility Research Demo ===\n\n";
    std::cout << "Fractional SDE terminal state: " << fractionalPath.states.back() << "\n\n";

    std::cout << "Target rough Bergomi parameters\n";
    std::cout << "  xi0=" << params.xi0
              << " eta=" << params.eta
              << " rho=" << params.rho
              << " H=" << params.hurst << "\n\n";

    std::cout << "Calibrated rough Bergomi parameters\n";
    std::cout << "  xi0=" << fit.parameters.xi0
              << " eta=" << fit.parameters.eta
              << " rho=" << fit.parameters.rho
              << " H=" << fit.parameters.hurst << "\n";
    std::cout << "  objective=SVI-smoothed surface RMSE\n";
    std::cout << "  surface RMSE=" << fit.rmse << "\n\n";

    std::cout << "Rough Bergomi pricing: MC vs QMC vs RQMC\n";
    std::cout << "  MC   price=" << mcPrice.price
              << " stdErr=" << mcPrice.stdErr << "\n";
    std::cout << "  QMC  price=" << qmcPrice.price
              << " stdErr=n/a\n";
    std::cout << "  RQMC price=" << rqmcPrice.price
              << " stdErr=" << rqmcPrice.stdErr << "\n\n";

    std::cout << "Correlated drivers check\n";
    std::cout << "  target rho=" << corr.targetRho
              << " empirical corr(dW_vol, dB_price)=" << corr.empiricalRho
              << " using " << corr.samples << " increments\n\n";

    std::cout << "Smile diagnostics at T=1.0\n";
    for (std::size_t i = 0; i < smile.size(); ++i) {
        std::cout << "  K=" << smile[i].strike << " IV=" << smile[i].impliedVol << "\n";
    }
    std::cout << "  Left-skew metric (low-strike IV - high-strike IV): " << leftSkew << "\n";
    std::cout << "  Equity-style left skew: " << (leftSkew > 0.0 ? "yes" : "no") << "\n\n";

    std::cout << "ATM term structure (K=S0)\n";
    for (std::size_t i = 0; i < term.size(); ++i) {
        std::cout << "  T=" << term[i].expiry << " IV=" << term[i].impliedVol << "\n";
    }
    std::cout << "\n";

    std::cout << "SVI comparison vs rough surface\n";
    std::cout << "  T=" << sviExpiry
              << " slice RMSE=" << sviRmse
              << " maxAbsErr=" << sviMaxErr << "\n";
    std::cout << "\n";

    std::cout << "Monte Carlo variance analysis\n";
    for (std::size_t i = 0; i < variance.size(); ++i) {
        std::cout << "  paths=" << variance[i].paths
                  << " RMSE=" << variance[i].rmse
                  << " stddev=" << variance[i].rmseStdDev << "\n";
    }
    std::cout << "  log-log slope (RMSE vs paths): " << logLogSlope(variance) << "\n\n";

    std::cout << "Diagnostics exported to " << OUTPUT_DIR << "\n";
    std::cout << "Run: python " << OUTPUT_DIR << "/plot_roughvol_diagnostics.py\n";

    return 0;
}
