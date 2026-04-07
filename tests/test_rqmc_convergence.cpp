#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

#include "../src/core/BrownianBridge.hpp"
#include "../src/core/Random.hpp"
#include "../src/core/SobolSequence.hpp"
#include "../src/stochasticProcess/GeometricBrownianMotion.hpp"

namespace {

constexpr double S0 = 100.0;
constexpr double K = 100.0;
constexpr double R = 0.05;
constexpr double SIGMA = 0.20;
constexpr double T = 1.0;
constexpr int RQMC_REPLICATES = 16;
constexpr int RANDOMIZED_OUTER_RUNS = 2;
constexpr std::uint64_t BASE_SEED = 0x1234ABCDEFULL;

const std::string DIR_FILE = "docs/new-joe-kuo-6.21201.txt";
const std::string OUTPUT_DIR = "tests/output";
const std::string CSV_FILE = OUTPUT_DIR + "/rqmc_convergence.csv";
const std::string PLOT_SCRIPT = OUTPUT_DIR + "/plot_rqmc_convergence.py";

const std::vector<int> PATH_COUNTS = {256, 512, 1024, 2048};
const std::vector<int> DIMENSIONS = {8, 32, 128};

struct Config {
    std::string method;
    std::string scramble;
    bool useBrownianBridge;
    bool useScrambling;
    bool includeInPlot;
    std::string label() const {
        std::ostringstream oss;
        oss << method;
        if (method == "RQMC") {
            oss << " (" << scramble << ")";
        }
        if (method == "MC") {
            oss << " (pseudo-random)";
        }
        return oss.str();
    }
};

struct Estimate {
    double price = 0.0;
    double error = 0.0;
    double stdErr = 0.0;
    double elapsedMs = 0.0;
};

struct DimensionOrderingMetrics {
    std::vector<double> orderedVarianceShare;
    int effectiveDimension50 = 0;
    int effectiveDimension90 = 0;
    int effectiveDimension99 = 0;
};

double rootMeanSquareError(const std::vector<double>& values, double target) {
    double sq = 0.0;
    for (double value : values) {
        const double d = value - target;
        sq += d * d;
    }
    return std::sqrt(sq / values.size());
}

double normalInverse(double u) {
    if (u <= 0.0) u = 1e-15;
    if (u >= 1.0) u = 1.0 - 1e-15;

    static const double a[] = {
        -3.969683028665376e+01,  2.209460984245205e+02,
        -2.759285104469687e+02,  1.383577518672690e+02,
        -3.066479806614716e+01,  2.506628277459239e+00
    };
    static const double b[] = {
        -5.447609879822406e+01,  1.615858368580409e+02,
        -1.556989798598866e+02,  6.680131188771972e+01,
        -1.328068155288572e+01
    };
    static const double c[] = {
        -7.784894002430293e-03, -3.223964580411365e-01,
        -2.400758277161838e+00, -2.549732539343734e+00,
         4.374664141464968e+00,  2.938163982698783e+00
    };
    static const double d[] = {
         7.784695709041462e-03,  3.224671290700398e-01,
         2.445134137142996e+00,  3.754408661907416e+00
    };

    const double pLow = 0.02425;
    const double pHigh = 1.0 - pLow;

    double q = 0.0;
    double r = 0.0;

    if (u < pLow) {
        q = std::sqrt(-2.0 * std::log(u));
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
    }

    if (u <= pHigh) {
        q = u - 0.5;
        r = q * q;
        return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
               (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
    }

    q = std::sqrt(-2.0 * std::log(1.0 - u));
    return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
            ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
}

std::uint64_t splitmix64(std::uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

double standardNormalCdf(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

double geometricAsianExactPrice(int steps) {
    const double dt = T / steps;
    const double avgTime = T * (steps + 1.0) / (2.0 * steps);
    const double variance = SIGMA * SIGMA * T * (steps + 1.0) * (2.0 * steps + 1.0)
                          / (6.0 * steps * steps);
    const double meanLog = std::log(S0) + (R - 0.5 * SIGMA * SIGMA) * avgTime;
    const double vol = std::sqrt(variance);

    const double d1 = (meanLog + variance - std::log(K)) / vol;
    const double d2 = d1 - vol;

    return std::exp(-R * T) *
           (std::exp(meanLog + 0.5 * variance) * standardNormalCdf(d1)
            - K * standardNormalCdf(d2));
}

double continuousGeometricAsianExactPrice() {
    const double avgTime = 0.5 * T;
    const double variance = SIGMA * SIGMA * T / 3.0;
    const double meanLog = std::log(S0) + (R - 0.5 * SIGMA * SIGMA) * avgTime;
    const double vol = std::sqrt(variance);

    const double d1 = (meanLog + variance - std::log(K)) / vol;
    const double d2 = d1 - vol;

    return std::exp(-R * T) *
           (std::exp(meanLog + 0.5 * variance) * standardNormalCdf(d1)
            - K * standardNormalCdf(d2));
}

double geometricAsianPayoff(
    int steps,
    const std::vector<double>& z,
    bool useBrownianBridge)
{
    const double dt = T / steps;
    const double sqrtDt = std::sqrt(dt);

    std::vector<double> pathNormals = z;
    if (useBrownianBridge) {
        auto dW = BrownianBridge::build(z, dt);
        for (int i = 0; i < steps; ++i) {
            pathNormals[i] = dW[i] / sqrtDt;
        }
    }

    double x = S0;
    double sumLog = 0.0;
    for (int i = 0; i < steps; ++i) {
        x *= std::exp((R - 0.5 * SIGMA * SIGMA) * dt + SIGMA * sqrtDt * pathNormals[i]);
        sumLog += std::log(std::max(x, 1e-12));
    }

    const double geometricAverage = std::exp(sumLog / steps);
    return std::exp(-R * T) * std::max(geometricAverage - K, 0.0);
}

double slopeFromLogLog(const std::vector<int>& x, const std::vector<double>& y) {
    std::vector<double> xs;
    std::vector<double> ys;

    for (std::size_t i = 0; i < x.size(); ++i) {
        if (y[i] > 0.0) {
            xs.push_back(std::log(static_cast<double>(x[i])));
            ys.push_back(std::log(y[i]));
        }
    }

    if (xs.size() < 2) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double meanX = std::accumulate(xs.begin(), xs.end(), 0.0) / xs.size();
    const double meanY = std::accumulate(ys.begin(), ys.end(), 0.0) / ys.size();

    double numerator = 0.0;
    double denominator = 0.0;
    for (std::size_t i = 0; i < xs.size(); ++i) {
        numerator += (xs[i] - meanX) * (ys[i] - meanY);
        denominator += (xs[i] - meanX) * (xs[i] - meanX);
    }

    return numerator / denominator;
}

std::vector<std::vector<double>> brownianBridgeMatrix(int steps, double dt) {
    std::vector<std::vector<double>> matrix(steps, std::vector<double>(steps, 0.0));
    for (int k = 0; k < steps; ++k) {
        std::vector<double> basis(steps, 0.0);
        basis[k] = 1.0;
        const auto dW = BrownianBridge::build(basis, dt);
        const double invSqrtDt = 1.0 / std::sqrt(dt);
        for (int j = 0; j < steps; ++j)
            matrix[j][k] = dW[j] * invSqrtDt;
    }
    return matrix;
}

DimensionOrderingMetrics analyzeDimensionOrdering(int steps, bool useBrownianBridge) {
    const double dt = T / steps;
    std::vector<double> incrementWeights(steps, 0.0);
    for (int j = 0; j < steps; ++j) {
        incrementWeights[j] = static_cast<double>(steps - j) / static_cast<double>(steps);
    }

    std::vector<double> sobolWeights(steps, 0.0);
    if (useBrownianBridge) {
        const auto transform = brownianBridgeMatrix(steps, dt);
        for (int k = 0; k < steps; ++k) {
            double coeff = 0.0;
            for (int j = 0; j < steps; ++j)
                coeff += incrementWeights[j] * transform[j][k];
            sobolWeights[k] = coeff * coeff;
        }
    } else {
        for (int k = 0; k < steps; ++k)
            sobolWeights[k] = incrementWeights[k] * incrementWeights[k];
    }

    const double total = std::accumulate(sobolWeights.begin(), sobolWeights.end(), 0.0);
    DimensionOrderingMetrics metrics;
    metrics.orderedVarianceShare.resize(steps, 0.0);
    double cumulative = 0.0;
    for (int k = 0; k < steps; ++k) {
        cumulative += sobolWeights[k];
        metrics.orderedVarianceShare[k] = (total > 0.0) ? cumulative / total : 0.0;
        if (metrics.effectiveDimension50 == 0 && metrics.orderedVarianceShare[k] >= 0.50)
            metrics.effectiveDimension50 = k + 1;
        if (metrics.effectiveDimension90 == 0 && metrics.orderedVarianceShare[k] >= 0.90)
            metrics.effectiveDimension90 = k + 1;
        if (metrics.effectiveDimension99 == 0 && metrics.orderedVarianceShare[k] >= 0.99)
            metrics.effectiveDimension99 = k + 1;
    }

    if (metrics.effectiveDimension50 == 0) metrics.effectiveDimension50 = steps;
    if (metrics.effectiveDimension90 == 0) metrics.effectiveDimension90 = steps;
    if (metrics.effectiveDimension99 == 0) metrics.effectiveDimension99 = steps;
    return metrics;
}

Estimate runMcEstimate(int steps, int totalPaths, bool useBrownianBridge, std::uint64_t seed, double exactPrice) {
    Random rng(static_cast<unsigned int>(seed & 0xFFFFFFFFULL));
    std::vector<double> payoffs;
    payoffs.reserve(totalPaths);
    std::vector<double> z(steps);

    const auto start = std::chrono::high_resolution_clock::now();
    for (int path = 0; path < totalPaths; ++path) {
        for (int j = 0; j < steps; ++j) {
            z[j] = rng.normal();
        }
        payoffs.push_back(geometricAsianPayoff(steps, z, useBrownianBridge));
    }
    const auto end = std::chrono::high_resolution_clock::now();

    const double mean = std::accumulate(payoffs.begin(), payoffs.end(), 0.0) / payoffs.size();
    double sq = 0.0;
    for (double p : payoffs) {
        const double d = p - mean;
        sq += d * d;
    }

    Estimate estimate;
    estimate.price = mean;
    estimate.error = std::abs(mean - exactPrice);
    estimate.stdErr = std::sqrt(sq / std::max(1, totalPaths - 1)) / std::sqrt(static_cast<double>(totalPaths));
    estimate.elapsedMs = std::chrono::duration<double, std::milli>(end - start).count();
    return estimate;
}

Estimate runQmcEstimate(int steps, int totalPaths, bool useBrownianBridge, double exactPrice) {
    SobolSequence sobol(steps, DIR_FILE);
    std::vector<double> z(steps);

    const auto start = std::chrono::high_resolution_clock::now();
    double sum = 0.0;
    for (int path = 0; path < totalPaths; ++path) {
        const auto u = sobol.next();
        for (int j = 0; j < steps; ++j) {
            z[j] = normalInverse(u[j]);
        }
        sum += geometricAsianPayoff(steps, z, useBrownianBridge);
    }
    const auto end = std::chrono::high_resolution_clock::now();

    Estimate estimate;
    estimate.price = sum / totalPaths;
    estimate.error = std::abs(estimate.price - exactPrice);
    estimate.stdErr = std::numeric_limits<double>::quiet_NaN();
    estimate.elapsedMs = std::chrono::duration<double, std::milli>(end - start).count();
    return estimate;
}

Estimate runRqmcEstimate(
    int steps,
    int totalPaths,
    bool useBrownianBridge,
    bool useScrambling,
    std::uint64_t seed,
    double exactPrice)
{
    if (totalPaths % RQMC_REPLICATES != 0) {
        throw std::runtime_error("RQMC totalPaths must be divisible by the replicate count");
    }

    const int pointsPerReplicate = totalPaths / RQMC_REPLICATES;

    std::vector<double> replicateMeans(RQMC_REPLICATES, 0.0);
    std::vector<double> shift(steps, 0.0);
    std::vector<double> z(steps);

    const auto start = std::chrono::high_resolution_clock::now();

    for (int rep = 0; rep < RQMC_REPLICATES; ++rep) {
        const std::uint64_t repSeed = seed + static_cast<std::uint64_t>(rep) * 0x9E3779B97F4A7C15ULL;
        SobolSequence sobol = useScrambling
            ? SobolSequence(steps, DIR_FILE, repSeed, 12)
            : SobolSequence(steps, DIR_FILE);

        if (!useScrambling) {
            std::uint64_t shiftSeed = splitmix64(repSeed);
            for (int j = 0; j < steps; ++j) {
                shiftSeed = splitmix64(shiftSeed);
                constexpr double inv2_53 = 1.0 / 9007199254740992.0;
                shift[j] = static_cast<double>((shiftSeed >> 11) & ((1ULL << 53) - 1)) * inv2_53;
            }
        }

        double payoffSum = 0.0;
        for (int n = 0; n < pointsPerReplicate; ++n) {
            auto u = sobol.next();
            for (int j = 0; j < steps; ++j) {
                double uj = u[j];
                if (!useScrambling) {
                    uj += shift[j];
                    if (uj >= 1.0) uj -= 1.0;
                }
                z[j] = normalInverse(uj);
            }
            payoffSum += geometricAsianPayoff(steps, z, useBrownianBridge);
        }
        replicateMeans[rep] = payoffSum / pointsPerReplicate;
    }

    const auto end = std::chrono::high_resolution_clock::now();

    const double mean = std::accumulate(replicateMeans.begin(), replicateMeans.end(), 0.0) / replicateMeans.size();
    double sq = 0.0;
    for (double value : replicateMeans) {
        const double d = value - mean;
        sq += d * d;
    }

    Estimate estimate;
    estimate.price = mean;
    estimate.error = std::abs(mean - exactPrice);
    estimate.stdErr = std::sqrt(sq / std::max(1, RQMC_REPLICATES - 1)) / std::sqrt(static_cast<double>(RQMC_REPLICATES));
    estimate.elapsedMs = std::chrono::duration<double, std::milli>(end - start).count();
    return estimate;
}

Estimate evaluateConfig(const Config& config, int steps, int totalPaths, double exactPrice) {
    if (config.method == "MC") {
        std::vector<double> prices;
        std::vector<double> innerStdErrs;
        double totalElapsed = 0.0;
        for (int run = 0; run < RANDOMIZED_OUTER_RUNS; ++run) {
            Estimate estimate = runMcEstimate(
                steps,
                totalPaths,
                false,
                BASE_SEED + static_cast<std::uint64_t>(steps * 100000 + totalPaths * 31 + run * 7919),
                exactPrice);
            prices.push_back(estimate.price);
            innerStdErrs.push_back(estimate.stdErr);
            totalElapsed += estimate.elapsedMs;
        }
        const double meanPrice = std::accumulate(prices.begin(), prices.end(), 0.0) / prices.size();
        const double meanStdErr = std::accumulate(innerStdErrs.begin(), innerStdErrs.end(), 0.0) / innerStdErrs.size();
        const double meanElapsed = totalElapsed / RANDOMIZED_OUTER_RUNS;
        const double rmsError = rootMeanSquareError(prices, exactPrice);
        return {meanPrice, rmsError, meanStdErr, meanElapsed};
    }

    if (config.method == "QMC") {
        return runQmcEstimate(steps, totalPaths, config.useBrownianBridge, exactPrice);
    }

    std::vector<double> prices;
    std::vector<double> innerStdErrs;
    double totalElapsed = 0.0;
    for (int run = 0; run < RANDOMIZED_OUTER_RUNS; ++run) {
        Estimate estimate = runRqmcEstimate(
            steps,
            totalPaths,
            config.useBrownianBridge,
            config.useScrambling,
            BASE_SEED + static_cast<std::uint64_t>(steps * 100000 + totalPaths * 53
                     + run * 1297 + (config.useScrambling ? 17 : 3)),
            exactPrice);
        prices.push_back(estimate.price);
        innerStdErrs.push_back(estimate.stdErr);
        totalElapsed += estimate.elapsedMs;
    }

    const double meanPrice = std::accumulate(prices.begin(), prices.end(), 0.0) / prices.size();
    const double meanStdErr = std::accumulate(innerStdErrs.begin(), innerStdErrs.end(), 0.0) / innerStdErrs.size();
    const double meanElapsed = totalElapsed / RANDOMIZED_OUTER_RUNS;
    const double rmsError = rootMeanSquareError(prices, exactPrice);
    return {meanPrice, rmsError, meanStdErr, meanElapsed};
}

std::string formatMetric(double value) {
    if (!std::isfinite(value)) {
        return "n/a";
    }

    std::ostringstream oss;
    oss << std::scientific << std::setprecision(4) << value;
    return oss.str();
}

std::string errorKindForMethod(const std::string& method) {
    return method == "QMC" ? "absolute_bias" : "rms_error";
}

std::string plotLabelForConfig(const Config& config) {
    if (config.method == "RQMC") {
        return "RQMC (" + config.scramble + ")";
    }
    return config.method;
}

void writePlotScript() {
    std::ofstream py(PLOT_SCRIPT);
    py << "import csv\n";
    py << "import math\n";
    py << "from pathlib import Path\n\n";
    py << "csv_path = Path(r'" << CSV_FILE << "')\n";
    py << R"PY(
out_svg = csv_path.with_suffix('.svg')
decomp_svg = csv_path.with_name('rqmc_error_decomposition.svg')
effdim_svg = csv_path.with_name('rqmc_effdim_vs_slope.svg')

rows = []
with csv_path.open(newline='') as f:
    reader = csv.DictReader(f)
    for row in reader:
        if not row or row.get('include_in_plot') in (None, ''):
            continue
        required = ['dimension', 'bridge', 'N', 'error', 'slope', 'sampling_error', 'discretization_bias', 'effective_dim_90']
        if any(row.get(key) in (None, '') for key in required):
            continue
        if int(row['include_in_plot']) != 1:
            continue
        row['dimension'] = int(row['dimension'])
        row['bridge'] = int(row['bridge'])
        row['N'] = int(row['N'])
        row['error'] = float(row['error'])
        row['slope'] = float(row['slope'])
        row['sampling_error'] = float(row['sampling_error'])
        row['discretization_bias'] = float(row['discretization_bias'])
        row['effective_dim_90'] = int(row['effective_dim_90'])
        rows.append(row)

dimensions = sorted({row['dimension'] for row in rows})
plot_order = ['MC', 'QMC', 'RQMC (Digital Shift)', 'RQMC (Owen Scramble)']
colors = {
    'MC': '#222222',
    'QMC': '#1f77b4',
    'RQMC (Digital Shift)': '#2ca02c',
    'RQMC (Owen Scramble)': '#d62728',
}

panel_w = 620
panel_h = 280
left_pad = 80
right_pad = 20
top_pad = 45
bottom_pad = 55
plot_w = panel_w - left_pad - right_pad
plot_h = panel_h - top_pad - bottom_pad
outer_pad = 20
svg_w = outer_pad * 3 + panel_w * 2
svg_h = outer_pad * (len(dimensions) + 1) + panel_h * len(dimensions)

def esc(text):
    return (str(text)
        .replace('&', '&amp;')
        .replace('<', '&lt;')
        .replace('>', '&gt;'))

def log_map(value, lo, hi, start, span):
    return start + (math.log10(value) - math.log10(lo)) / (math.log10(hi) - math.log10(lo)) * span

def grouped(panel_rows):
    by_label = {}
    for label in plot_order:
        series = [row for row in panel_rows if row['plot_label'] == label]
        if series:
            by_label[label] = sorted(series, key=lambda row: row['N'])
    return by_label

parts = [
    f'<svg xmlns="http://www.w3.org/2000/svg" width="{svg_w}" height="{svg_h}" viewBox="0 0 {svg_w} {svg_h}">',
    '<style>',
    'text { font-family: Arial, sans-serif; fill: #111; }',
    '.title { font-size: 20px; font-weight: bold; }',
    '.subtitle { font-size: 12px; fill: #444; }',
    '.panel-title { font-size: 14px; font-weight: bold; }',
    '.axis { font-size: 11px; fill: #444; }',
    '.legend { font-size: 11px; }',
    '.grid { stroke: #cccccc; stroke-width: 1; opacity: 0.45; }',
    '.frame { stroke: #444; stroke-width: 1; fill: none; }',
    '.guide { stroke-width: 1.2; fill: none; opacity: 0.85; }',
    '.series { stroke-width: 2.2; fill: none; }',
    '</style>',
]

parts.append(f'<text x="{svg_w / 2}" y="26" text-anchor="middle" class="title">MC vs QMC vs RQMC convergence</text>')
parts.append(f'<text x="{svg_w / 2}" y="44" text-anchor="middle" class="subtitle">Geometric-Asian call under GBM, log-log error vs total paths N</text>')

for row_idx, dim in enumerate(dimensions):
    for col_idx, bridge in enumerate([0, 1]):
        x0 = outer_pad + col_idx * (panel_w + outer_pad)
        y0 = outer_pad + row_idx * (panel_h + outer_pad)
        plot_x0 = x0 + left_pad
        plot_y0 = y0 + top_pad

        panel_rows = [row for row in rows if row['dimension'] == dim and row['bridge'] == bridge]
        series_map = grouped(panel_rows)
        all_x = sorted({row['N'] for row in panel_rows})
        all_y = [row['error'] for row in panel_rows if row['error'] > 0.0]
        if not all_x or not all_y:
            continue

        y_min = min(all_y)
        y_max = max(all_y)
        y_min = 10 ** math.floor(math.log10(y_min))
        y_max = 10 ** math.ceil(math.log10(y_max))
        if y_min == y_max:
            y_min /= 10.0
            y_max *= 10.0

        parts.append(f'<rect x="{x0}" y="{y0}" width="{panel_w}" height="{panel_h}" fill="#ffffff" stroke="#dddddd"/>')
        bridge_label = 'ON' if bridge else 'OFF'
        parts.append(f'<text x="{x0 + panel_w / 2}" y="{y0 + 20}" text-anchor="middle" class="panel-title">d = {dim}, Brownian bridge {bridge_label}</text>')

        for decade in range(int(math.log10(y_min)), int(math.log10(y_max)) + 1):
            tick_val = 10 ** decade
            yy = plot_y0 + plot_h - log_map(tick_val, y_min, y_max, 0, plot_h)
            parts.append(f'<line x1="{plot_x0}" y1="{yy:.2f}" x2="{plot_x0 + plot_w}" y2="{yy:.2f}" class="grid"/>')
            parts.append(f'<text x="{plot_x0 - 8}" y="{yy + 4:.2f}" text-anchor="end" class="axis">1e{decade}</text>')

        for xv in all_x:
            xx = log_map(xv, min(all_x), max(all_x), plot_x0, plot_w)
            parts.append(f'<line x1="{xx:.2f}" y1="{plot_y0}" x2="{xx:.2f}" y2="{plot_y0 + plot_h}" class="grid"/>')
            parts.append(f'<text x="{xx:.2f}" y="{plot_y0 + plot_h + 18}" text-anchor="middle" class="axis">{xv}</text>')

        parts.append(f'<rect x="{plot_x0}" y="{plot_y0}" width="{plot_w}" height="{plot_h}" class="frame"/>')
        parts.append(f'<text x="{x0 + panel_w / 2}" y="{y0 + panel_h - 12}" text-anchor="middle" class="axis">N (total paths)</text>')
        parts.append(f'<text x="{x0 + 18}" y="{y0 + panel_h / 2}" text-anchor="middle" transform="rotate(-90 {x0 + 18} {y0 + panel_h / 2})" class="axis">Absolute error</text>')

        if 'MC' in series_map:
            mc = series_map['MC']
            x_ref = mc[0]['N']
            y_ref = mc[0]['error']
            points = []
            for xv in all_x:
                yv = y_ref * (xv / x_ref) ** (-0.5)
                xx = log_map(xv, min(all_x), max(all_x), plot_x0, plot_w)
                yy = plot_y0 + plot_h - log_map(yv, y_min, y_max, 0, plot_h)
                points.append(f'{xx:.2f},{yy:.2f}')
            parts.append(f'<polyline points="{" ".join(points)}" class="guide" stroke="#666666" stroke-dasharray="6,4"/>')

        if 'QMC' in series_map:
            qmc = series_map['QMC']
            x_ref = qmc[0]['N']
            y_ref = qmc[0]['error']
            points = []
            for xv in all_x:
                yv = y_ref * (xv / x_ref) ** (-1.0)
                xx = log_map(xv, min(all_x), max(all_x), plot_x0, plot_w)
                yy = plot_y0 + plot_h - log_map(yv, y_min, y_max, 0, plot_h)
                points.append(f'{xx:.2f},{yy:.2f}')
            parts.append(f'<polyline points="{" ".join(points)}" class="guide" stroke="#8888cc" stroke-dasharray="2,4"/>')

        for label in plot_order:
            if label not in series_map:
                continue
            series = series_map[label]
            pts = []
            for item in series:
                xx = log_map(item['N'], min(all_x), max(all_x), plot_x0, plot_w)
                yy = plot_y0 + plot_h - log_map(item['error'], y_min, y_max, 0, plot_h)
                pts.append((xx, yy))
            point_str = ' '.join(f'{xx:.2f},{yy:.2f}' for xx, yy in pts)
            parts.append(f'<polyline points="{point_str}" class="series" stroke="{colors[label]}"/>')
            for xx, yy in pts:
                parts.append(f'<circle cx="{xx:.2f}" cy="{yy:.2f}" r="3.2" fill="{colors[label]}"/>')

        legend_x = plot_x0 + 8
        legend_y = plot_y0 + 16
        legend_items = []
        for label in plot_order:
            if label in series_map:
                slope = series_map[label][0]['slope']
                legend_items.append((label, f'{label} (slope={slope:.2f})'))
        legend_items.append(('guide_mc', 'MC guide N^-1/2'))
        legend_items.append(('guide_qmc', 'Ideal QMC guide N^-1'))

        for idx, item in enumerate(legend_items):
            key, text = item
            yy = legend_y + idx * 16
            if key == 'guide_mc':
                parts.append(f'<line x1="{legend_x}" y1="{yy}" x2="{legend_x + 18}" y2="{yy}" stroke="#666666" stroke-width="1.2" stroke-dasharray="6,4"/>')
            elif key == 'guide_qmc':
                parts.append(f'<line x1="{legend_x}" y1="{yy}" x2="{legend_x + 18}" y2="{yy}" stroke="#8888cc" stroke-width="1.2" stroke-dasharray="2,4"/>')
            else:
                parts.append(f'<line x1="{legend_x}" y1="{yy}" x2="{legend_x + 18}" y2="{yy}" stroke="{colors[key]}" stroke-width="2.2"/>')
                parts.append(f'<circle cx="{legend_x + 9}" cy="{yy}" r="3.0" fill="{colors[key]}"/>')
            parts.append(f'<text x="{legend_x + 24}" y="{yy + 4}" class="legend">{esc(text)}</text>')

parts.append('</svg>')
out_svg.write_text('\n'.join(parts), encoding='utf-8')
print(f'Wrote {out_svg}')

decomp_rows = [row for row in rows if row['N'] == max(r['N'] for r in rows)]
bar_w = 180
bar_h = 240
svg_w2 = outer_pad * 3 + bar_w * 2
svg_h2 = outer_pad * (len(dimensions) + 1) + bar_h * len(dimensions)
parts2 = [
    f'<svg xmlns="http://www.w3.org/2000/svg" width="{svg_w2}" height="{svg_h2}" viewBox="0 0 {svg_w2} {svg_h2}">',
    '<style>',
    'text { font-family: Arial, sans-serif; fill: #111; }',
    '.title { font-size: 20px; font-weight: bold; }',
    '.axis { font-size: 11px; fill: #444; }',
    '.label { font-size: 11px; }',
    '.frame { stroke: #444; stroke-width: 1; fill: none; }',
    '</style>',
]
parts2.append(f'<text x="{svg_w2 / 2}" y="26" text-anchor="middle" class="title">QMC Error Decomposition</text>')
max_err = max(max(r['sampling_error'], r['discretization_bias']) for r in decomp_rows)
max_err = 1.1 * max_err if max_err > 0 else 1.0
for row_idx, dim in enumerate(dimensions):
    for col_idx, bridge in enumerate([0, 1]):
        panel = [r for r in decomp_rows if r['dimension'] == dim and r['bridge'] == bridge and r['plot_label'] == 'QMC']
        if not panel:
            continue
        item = panel[0]
        x0 = outer_pad + col_idx * (bar_w + outer_pad)
        y0 = outer_pad + row_idx * (bar_h + outer_pad)
        parts2.append(f'<rect x="{x0}" y="{y0}" width="{bar_w}" height="{bar_h}" fill="#ffffff" stroke="#dddddd"/>')
        parts2.append(f'<text x="{x0 + bar_w / 2}" y="{y0 + 18}" text-anchor="middle" class="label">d={dim}, bridge {"ON" if bridge else "OFF"}</text>')
        plot_y0 = y0 + 30
        plot_h2 = 140
        base_y = plot_y0 + plot_h2
        sampling_h = plot_h2 * item['sampling_error'] / max_err
        bias_h = plot_h2 * item['discretization_bias'] / max_err
        parts2.append(f'<rect x="{x0 + 30}" y="{base_y - sampling_h:.2f}" width="36" height="{sampling_h:.2f}" fill="#1f77b4"/>')
        parts2.append(f'<rect x="{x0 + 88}" y="{base_y - bias_h:.2f}" width="36" height="{bias_h:.2f}" fill="#d62728"/>')
        parts2.append(f'<line x1="{x0 + 18}" y1="{base_y}" x2="{x0 + bar_w - 18}" y2="{base_y}" class="frame"/>')
        parts2.append(f'<text x="{x0 + 48}" y="{base_y + 16}" text-anchor="middle" class="axis">sampling</text>')
        parts2.append(f'<text x="{x0 + 106}" y="{base_y + 16}" text-anchor="middle" class="axis">bias</text>')
        parts2.append(f'<text x="{x0 + bar_w / 2}" y="{y0 + 195}" text-anchor="middle" class="label">effective dim 90% = {item["effective_dim_90"]}</text>')
        parts2.append(f'<text x="{x0 + bar_w / 2}" y="{y0 + 212}" text-anchor="middle" class="label">ordering impact shown by early-dim concentration</text>')
parts2.append('</svg>')
decomp_svg.write_text('\n'.join(parts2), encoding='utf-8')
print(f'Wrote {decomp_svg}')

slope_rows = []
seen = set()
for row in rows:
    key = (row['dimension'], row['bridge'], row['plot_label'])
    if key in seen:
        continue
    seen.add(key)
    slope_rows.append(row)

qmc_like = [row for row in slope_rows if row['plot_label'] != 'MC']
if qmc_like:
    min_x = min(row['effective_dim_90'] for row in qmc_like)
    max_x = max(row['effective_dim_90'] for row in qmc_like)
    if min_x == max_x:
        min_x -= 1
        max_x += 1

    min_y = min(min(row['slope'] for row in qmc_like), -1.0)
    max_y = max(max(row['slope'] for row in qmc_like), -0.5)
    y_pad = 0.08 * max(0.2, max_y - min_y)
    min_y -= y_pad
    max_y += y_pad

    svg_w3 = 760
    svg_h3 = 520
    left3 = 90
    right3 = 30
    top3 = 55
    bottom3 = 70
    plot_w3 = svg_w3 - left3 - right3
    plot_h3 = svg_h3 - top3 - bottom3

    def x_map(x):
        return left3 + (x - min_x) / (max_x - min_x) * plot_w3

    def y_map(y):
        return top3 + plot_h3 - (y - min_y) / (max_y - min_y) * plot_h3

    marker_map = {
        (0, 'QMC'): 'circle',
        (1, 'QMC'): 'square',
        (0, 'RQMC (Digital Shift)'): 'triangle',
        (1, 'RQMC (Digital Shift)'): 'diamond',
        (0, 'RQMC (Owen Scramble)'): 'triangle_down',
        (1, 'RQMC (Owen Scramble)'): 'hex',
    }

    parts3 = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{svg_w3}" height="{svg_h3}" viewBox="0 0 {svg_w3} {svg_h3}">',
        '<style>',
        'text { font-family: Arial, sans-serif; fill: #111; }',
        '.title { font-size: 20px; font-weight: bold; }',
        '.axis { font-size: 12px; fill: #444; }',
        '.legend { font-size: 11px; }',
        '.grid { stroke: #cccccc; stroke-width: 1; opacity: 0.45; }',
        '.frame { stroke: #444; stroke-width: 1; fill: none; }',
        '.guide-qmc { stroke: #1f77b4; stroke-width: 1.5; stroke-dasharray: 4,4; opacity: 0.8; }',
        '.guide-mc { stroke: #666666; stroke-width: 1.5; stroke-dasharray: 6,4; opacity: 0.8; }',
        '</style>',
        f'<text x="{svg_w3 / 2}" y="28" text-anchor="middle" class="title">Effective Dimension vs QMC Convergence Slope</text>',
        f'<text x="{svg_w3 / 2}" y="46" text-anchor="middle" class="axis">As EffDim90 rises, QMC slope typically degrades from ideal -1 toward MC-like -0.5</text>',
    ]

    x_ticks = sorted(set(row['effective_dim_90'] for row in qmc_like))
    for x_tick in x_ticks:
        xx = x_map(x_tick)
        parts3.append(f'<line x1="{xx:.2f}" y1="{top3}" x2="{xx:.2f}" y2="{top3 + plot_h3}" class="grid"/>')
        parts3.append(f'<text x="{xx:.2f}" y="{top3 + plot_h3 + 20}" text-anchor="middle" class="axis">{x_tick}</text>')

    for y_tick in [-1.0, -0.9, -0.8, -0.7, -0.6, -0.5]:
        if y_tick < min_y or y_tick > max_y:
            continue
        yy = y_map(y_tick)
        parts3.append(f'<line x1="{left3}" y1="{yy:.2f}" x2="{left3 + plot_w3}" y2="{yy:.2f}" class="grid"/>')
        parts3.append(f'<text x="{left3 - 10}" y="{yy + 4:.2f}" text-anchor="end" class="axis">{y_tick:.1f}</text>')

    parts3.append(f'<rect x="{left3}" y="{top3}" width="{plot_w3}" height="{plot_h3}" class="frame"/>')
    parts3.append(f'<text x="{svg_w3 / 2}" y="{svg_h3 - 18}" text-anchor="middle" class="axis">Effective Dimension (EffDim90)</text>')
    parts3.append(f'<text x="24" y="{svg_h3 / 2}" text-anchor="middle" transform="rotate(-90 24 {svg_h3 / 2})" class="axis">Log-log convergence slope</text>')

    y_qmc = y_map(-1.0)
    y_mc = y_map(-0.5)
    parts3.append(f'<line x1="{left3}" y1="{y_qmc:.2f}" x2="{left3 + plot_w3}" y2="{y_qmc:.2f}" class="guide-qmc"/>')
    parts3.append(f'<line x1="{left3}" y1="{y_mc:.2f}" x2="{left3 + plot_w3}" y2="{y_mc:.2f}" class="guide-mc"/>')
    parts3.append(f'<text x="{left3 + plot_w3 - 8}" y="{y_qmc - 6:.2f}" text-anchor="end" class="legend">Ideal QMC slope -1</text>')
    parts3.append(f'<text x="{left3 + plot_w3 - 8}" y="{y_mc - 6:.2f}" text-anchor="end" class="legend">MC slope -0.5</text>')

    for row in qmc_like:
        xx = x_map(row['effective_dim_90'])
        yy = y_map(row['slope'])
        label = row['plot_label']
        color = colors.get(label, '#000000')
        marker = marker_map.get((row['bridge'], label), 'circle')
        if marker == 'circle':
            parts3.append(f'<circle cx="{xx:.2f}" cy="{yy:.2f}" r="5" fill="{color}"/>')
        elif marker == 'square':
            parts3.append(f'<rect x="{xx - 5:.2f}" y="{yy - 5:.2f}" width="10" height="10" fill="{color}"/>')
        elif marker == 'triangle':
            parts3.append(f'<polygon points="{xx:.2f},{yy - 6:.2f} {xx - 6:.2f},{yy + 5:.2f} {xx + 6:.2f},{yy + 5:.2f}" fill="{color}"/>')
        elif marker == 'diamond':
            parts3.append(f'<polygon points="{xx:.2f},{yy - 6:.2f} {xx - 6:.2f},{yy:.2f} {xx:.2f},{yy + 6:.2f} {xx + 6:.2f},{yy:.2f}" fill="{color}"/>')
        elif marker == 'triangle_down':
            parts3.append(f'<polygon points="{xx - 6:.2f},{yy - 5:.2f} {xx + 6:.2f},{yy - 5:.2f} {xx:.2f},{yy + 6:.2f}" fill="{color}"/>')
        else:
            parts3.append(f'<circle cx="{xx:.2f}" cy="{yy:.2f}" r="3.5" fill="{color}"/>')
            parts3.append(f'<circle cx="{xx:.2f}" cy="{yy:.2f}" r="6" fill="none" stroke="{color}" stroke-width="1.5"/>')
        bridge_label = 'bridge ON' if row['bridge'] == 1 else 'bridge OFF'
        parts3.append(f'<text x="{xx + 8:.2f}" y="{yy - 6:.2f}" class="legend">{esc(label)}, {bridge_label}, d={row["dimension"]}</text>')

    parts3.append('</svg>')
    effdim_svg.write_text('\n'.join(parts3), encoding='utf-8')
    print(f'Wrote {effdim_svg}')
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

    std::vector<Config> configs = {
        {"MC", "none", false, false, true},
        {"QMC", "none", false, false, true},
        {"QMC", "none", true, false, true},
        {"RQMC", "Digital Shift", false, false, true},
        {"RQMC", "Digital Shift", true, false, true},
        {"RQMC", "Owen Scramble", false, true, true},
        {"RQMC", "Owen Scramble", true, true, true},
    };

    std::ofstream csv(CSV_FILE);
    csv << "dimension,dt,N,method,scramble,bridge,price,error,sampling_error,discretization_bias,std_err,elapsed_ms,include_in_plot,plot_label,slope,error_kind,effective_dim_50,effective_dim_90,effective_dim_99\n";

    std::cout << "===============================================================\n";
    std::cout << "RQMC Convergence Study: MC vs QMC vs RQMC\n";
    std::cout << "Payoff: geometric-Asian call under GBM\n";
    std::cout << "Simulation: exact GBM step, so changing d does not add Euler bias\n";
    std::cout << "Theory guides: MC ~= O(N^{-1/2}), ideal QMC ~= O(N^{-1})\n";
    std::cout << "Brownian bridge: ON/OFF | Sobol scrambling: OFF(digital shift)/ON(Owen)\n";
    std::cout << "Error decomposition: total error = sampling error + discretization bias + ordering diagnostics\n";
    std::cout << "Sampling metric: RMS error for MC/RQMC, absolute bias for deterministic QMC against the discrete-time exact target\n";
    std::cout << "CSV output: " << CSV_FILE << "\n";
    std::cout << "Plot script: " << PLOT_SCRIPT << "\n";
    std::cout << "===============================================================\n\n";

    for (int steps : DIMENSIONS) {
        const double dt = T / steps;
        const double exactPrice = geometricAsianExactPrice(steps);
        const double continuousExact = continuousGeometricAsianExactPrice();
        const double discretizationBias = std::abs(exactPrice - continuousExact);

        std::cout << "Dimension d = " << steps
                  << " (dt=" << dt << ")"
                  << " | exact geometric-Asian price = "
                  << std::fixed << std::setprecision(8) << exactPrice << "\n";

        for (int bridgeState = 0; bridgeState <= 1; ++bridgeState) {
            std::cout << "  Brownian bridge " << (bridgeState ? "ON" : "OFF") << "\n";
            std::cout << std::left
                      << "  " << std::setw(24) << "Method"
                      << std::setw(12) << "Slope"
                      << std::setw(14) << "Samp@256"
                      << std::setw(14) << "Samp@2048"
                      << std::setw(14) << "DiscBias"
                      << std::setw(12) << "EffDim90"
                      << std::setw(14) << "StdErr@2048"
                      << "\n";
            std::cout << "  " << std::string(94, '-') << "\n";

            for (const Config& config : configs) {
                if (config.method != "MC" &&
                    config.useBrownianBridge != static_cast<bool>(bridgeState)) {
                    continue;
                }

                const bool bridgeForMetrics = (config.method == "MC") ? static_cast<bool>(bridgeState) : config.useBrownianBridge;
                const DimensionOrderingMetrics orderingMetrics = analyzeDimensionOrdering(steps, bridgeForMetrics);

                std::vector<double> errors;
                std::vector<Estimate> estimates;

                for (int totalPaths : PATH_COUNTS) {
                    Estimate estimate = evaluateConfig(config, steps, totalPaths, exactPrice);
                    estimates.push_back(estimate);
                    errors.push_back(estimate.error);
                }

                const double slope = slopeFromLogLog(PATH_COUNTS, errors);
                const std::string plotLabel = plotLabelForConfig(config);

                for (std::size_t i = 0; i < PATH_COUNTS.size(); ++i) {
                    csv << steps << ","
                        << dt << ","
                        << PATH_COUNTS[i] << ","
                        << config.method << ","
                        << config.scramble << ","
                        << (config.method == "MC" ? bridgeState : (config.useBrownianBridge ? 1 : 0)) << ","
                        << estimates[i].price << ","
                        << estimates[i].error << ","
                        << estimates[i].error << ","
                        << discretizationBias << ","
                        << (std::isfinite(estimates[i].stdErr) ? std::to_string(estimates[i].stdErr) : std::string("nan")) << ","
                        << estimates[i].elapsedMs << ","
                        << (config.includeInPlot ? 1 : 0) << ","
                        << plotLabel << ","
                        << slope << ","
                        << errorKindForMethod(config.method) << ","
                        << orderingMetrics.effectiveDimension50 << ","
                        << orderingMetrics.effectiveDimension90 << ","
                        << orderingMetrics.effectiveDimension99 << "\n";
                }

                std::string displayLabel = config.label();
                if (config.method == "MC") {
                    displayLabel += " [bridge-invariant]";
                }

                std::cout << "  " << std::setw(24) << displayLabel
                          << std::setw(12) << std::fixed << std::setprecision(3) << slope
                          << std::setw(14) << formatMetric(errors.front())
                          << std::setw(14) << formatMetric(errors.back())
                          << std::setw(14) << formatMetric(discretizationBias)
                          << std::setw(12) << orderingMetrics.effectiveDimension90
                          << std::setw(14) << formatMetric(estimates.back().stdErr)
                          << "\n";
            }
            std::cout << "\n";
        }
    }

    csv.close();

    std::cout << "Finished. Run the plot script with:\n";
    std::cout << "  python " << PLOT_SCRIPT << "\n";
    std::cout << "This produces convergence, decomposition, and effective-dimension plots in " << OUTPUT_DIR << ".\n";
    return 0;
}
