#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

#include "../src/integrators/EulerMaruyama.hpp"
#include "../src/payoffs/VanillaPayoffs.hpp"
#include "../src/pricing/EuropeanPricer.hpp"
#include "../src/simulators/ParallelMonteCarloSimulator.hpp"
#include "../src/stochasticProcess/GeometricBrownianMotion.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

const std::string OUTPUT_DIR = "tests/output";
const std::string CSV_FILE = OUTPUT_DIR + "/parallel_scaling.csv";
const std::string PLOT_SCRIPT = OUTPUT_DIR + "/plot_parallel_scaling.py";

struct ScalingPoint {
    std::string study;
    int threads;
    int paths;
    double meanTimeSec;
    double stddevTimeSec;
    double speedup;
    double efficiency;
    double price;
    double priceStdErr;
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

void ensureOutputDir() {
#ifdef _WIN32
    _mkdir(OUTPUT_DIR.c_str());
#else
    mkdir(OUTPUT_DIR.c_str(), 0755);
#endif
}

std::vector<int> threadGrid(int maxThreads) {
    std::vector<int> grid;
    for (int t = 1; t <= maxThreads; t *= 2)
        grid.push_back(t);
    if (grid.empty() || grid.back() != maxThreads)
        grid.push_back(maxThreads);
    grid.erase(std::unique(grid.begin(), grid.end()), grid.end());
    return grid;
}

ScalingPoint runPoint(
    const std::string& study,
    int threads,
    int paths,
    int repeats,
    GeometricBrownianMotion& gbm,
    double spot,
    double rate,
    double maturity,
    double dt,
    double strike)
{
    std::vector<double> times;
    times.reserve(repeats);
    double lastPrice = 0.0;
    double lastErr = 0.0;

    for (int rep = 0; rep < repeats; ++rep) {
        const auto t0 = std::chrono::high_resolution_clock::now();
        const auto terminal = ParallelMonteCarloSimulator<EulerMaruyama>::simulate(
            gbm, spot, maturity, dt, paths, threads);
        const auto t1 = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double> elapsed = t1 - t0;
        times.push_back(elapsed.count());

        const CallPayoff payoff(strike);
        const auto priced = EuropeanPricer::priceWithError(terminal, payoff, rate, maturity);
        lastPrice = priced.first;
        lastErr = priced.second;
    }

    ScalingPoint point;
    point.study = study;
    point.threads = threads;
    point.paths = paths;
    point.meanTimeSec = mean(times);
    point.stddevTimeSec = stdDev(times);
    point.speedup = 0.0;
    point.efficiency = 0.0;
    point.price = lastPrice;
    point.priceStdErr = lastErr;
    return point;
}

void writePlotScript() {
    std::ofstream py(PLOT_SCRIPT.c_str());
    py << "import csv\n";
    py << "from pathlib import Path\n\n";
    py << "csv_path = Path(r'" << CSV_FILE << "')\n";
    py << R"PY(
rows = list(csv.DictReader(csv_path.open(newline='')))
out_svg = csv_path.with_suffix('.svg')

studies = ['strong', 'weak']
colors = {'strong': '#1f77b4', 'weak': '#d62728'}
series = {s: {'threads': [], 'speedup': [], 'eff': [], 'time': [], 'paths': []} for s in studies}
for row in rows:
    s = row['study']
    if s not in series:
        continue
    series[s]['threads'].append(int(row['threads']))
    series[s]['speedup'].append(float(row['speedup']))
    series[s]['eff'].append(float(row['efficiency']))
    series[s]['time'].append(float(row['mean_time_sec']))
    series[s]['paths'].append(int(row['paths']))

def esc(text):
    return str(text).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')

def panel(parts, title, ylabel, key, x0, y0, w, h, guides=None):
    left, right, top, bottom = 64, 20, 38, 55
    pw = w - left - right
    ph = h - top - bottom
    px = x0 + left
    py = y0 + top
    xs = [x for s in studies for x in series[s]['threads']]
    ys = [y for s in studies for y in series[s][key]]
    if not xs or not ys:
        return
    if guides:
        ys = ys + guides
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    if xmax == xmin:
        xmax = xmin + 1
    if ymax == ymin:
        ymax = ymin + 1
    pad = 0.08 * (ymax - ymin)
    ymin -= pad
    ymax += pad
    parts.append(f'<rect x="{x0}" y="{y0}" width="{w}" height="{h}" fill="#fff" stroke="#dddddd"/>')
    parts.append(f'<text x="{x0 + w/2}" y="{y0 + 20}" text-anchor="middle" class="panel">{esc(title)}</text>')
    for i in range(5):
        xv = xmin + (xmax - xmin) * i / 4.0
        xx = px + (xv - xmin) / (xmax - xmin) * pw
        parts.append(f'<line x1="{xx:.2f}" y1="{py}" x2="{xx:.2f}" y2="{py + ph}" class="grid"/>')
        parts.append(f'<text x="{xx:.2f}" y="{py + ph + 18}" text-anchor="middle" class="axis">{xv:.0f}</text>')
    for i in range(5):
        yv = ymin + (ymax - ymin) * i / 4.0
        yy = py + ph - (yv - ymin) / (ymax - ymin) * ph
        parts.append(f'<line x1="{px}" y1="{yy:.2f}" x2="{px + pw}" y2="{yy:.2f}" class="grid"/>')
        parts.append(f'<text x="{px - 8}" y="{yy + 4:.2f}" text-anchor="end" class="axis">{yv:.3g}</text>')
    if guides:
        for gv in guides:
            yy = py + ph - (gv - ymin) / (ymax - ymin) * ph
            parts.append(f'<line x1="{px}" y1="{yy:.2f}" x2="{px + pw}" y2="{yy:.2f}" stroke="#888" stroke-dasharray="5,4" stroke-width="1.2"/>')
    parts.append(f'<rect x="{px}" y="{py}" width="{pw}" height="{ph}" class="frame"/>')
    for s in studies:
        pts = []
        for x, y in zip(series[s]['threads'], series[s][key]):
            xx = px + (x - xmin) / (xmax - xmin) * pw
            yy = py + ph - (y - ymin) / (ymax - ymin) * ph
            pts.append((xx, yy))
        if not pts:
            continue
        point_str = ' '.join(f'{xx:.2f},{yy:.2f}' for xx, yy in pts)
        parts.append(f'<polyline points="{point_str}" fill="none" stroke="{colors[s]}" stroke-width="2.2"/>')
        for xx, yy in pts:
            parts.append(f'<circle cx="{xx:.2f}" cy="{yy:.2f}" r="3" fill="{colors[s]}"/>')
    parts.append(f'<text x="{x0 + w/2}" y="{y0 + h - 14}" text-anchor="middle" class="axis">Threads / cores</text>')
    parts.append(f'<text x="{x0 + 18}" y="{y0 + h/2}" text-anchor="middle" transform="rotate(-90 {x0 + 18} {y0 + h/2})" class="axis">{esc(ylabel)}</text>')

parts = [
    '<svg xmlns="http://www.w3.org/2000/svg" width="1120" height="780" viewBox="0 0 1120 780">',
    '<style>',
    'text { font-family: Arial, sans-serif; fill: #111; }',
    '.title { font-size: 18px; font-weight: bold; }',
    '.panel { font-size: 13px; font-weight: bold; }',
    '.axis { font-size: 11px; fill: #444; }',
    '.legend { font-size: 11px; }',
    '.grid { stroke: #cccccc; stroke-width: 1; opacity: 0.45; }',
    '.frame { stroke: #444; stroke-width: 1; fill: none; }',
    '</style>',
    '<text x="560" y="24" text-anchor="middle" class="title">Parallel Scaling Study</text>',
]
panel(parts, 'Speedup vs Cores', 'Speedup', 'speedup', 20, 50, 520, 320, guides=[1.0])
panel(parts, 'Efficiency vs Cores', 'Efficiency', 'eff', 580, 50, 520, 320, guides=[1.0])
panel(parts, 'Runtime vs Cores', 'Mean wall time (s)', 'time', 300, 400, 520, 320)
lx = 90
ly = 86
labels = {'strong': 'Strong scaling', 'weak': 'Weak scaling'}
for i, s in enumerate(studies):
    yy = ly + i * 16
    parts.append(f'<line x1="{lx}" y1="{yy}" x2="{lx + 18}" y2="{yy}" stroke="{colors[s]}" stroke-width="2.2"/>')
    parts.append(f'<circle cx="{lx + 9}" cy="{yy}" r="3" fill="{colors[s]}"/>')
    parts.append(f'<text x="{lx + 24}" y="{yy + 4}" class="legend">{esc(labels[s])}</text>')
parts.append('</svg>')
out_svg.write_text('\n'.join(parts), encoding='utf-8')
print(f'Wrote {out_svg}')
)PY";
}

} // namespace

int main() {
    ensureOutputDir();
    writePlotScript();

#ifndef _OPENMP
    std::cerr << "This scaling study requires OpenMP. Rebuild with -fopenmp.\n";
    return 1;
#else
    GeometricBrownianMotion gbm(0.05, 0.2);
    const double spot = 100.0;
    const double rate = 0.05;
    const double maturity = 1.0;
    const double dt = 0.01;
    const double strike = 100.0;
    const int repeats = 3;
    const int strongPaths = 200000;
    const int weakPathsPerThread = 80000;

    const int maxThreads = std::max(1, omp_get_max_threads());
    const std::vector<int> threads = threadGrid(maxThreads);

    std::vector<ScalingPoint> points;
    points.reserve(threads.size() * 2);

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "=== Parallel Scaling Study ===\n";
    std::cout << "OpenMP max threads: " << maxThreads << "\n";

    for (int t : threads) {
        ScalingPoint point = runPoint(
            "strong", t, strongPaths, repeats,
            gbm, spot, rate, maturity, dt, strike);
        points.push_back(point);
        std::cout << "Strong scaling  threads=" << t
                  << " paths=" << strongPaths
                  << " time=" << point.meanTimeSec
                  << "s price=" << point.price
                  << " stderr=" << point.priceStdErr << "\n";
    }

    for (int t : threads) {
        ScalingPoint point = runPoint(
            "weak", t, weakPathsPerThread * t, repeats,
            gbm, spot, rate, maturity, dt, strike);
        points.push_back(point);
        std::cout << "Weak scaling    threads=" << t
                  << " paths=" << (weakPathsPerThread * t)
                  << " time=" << point.meanTimeSec
                  << "s price=" << point.price
                  << " stderr=" << point.priceStdErr << "\n";
    }

    double strongBase = 0.0;
    double weakBase = 0.0;
    for (const auto& p : points) {
        if (p.study == "strong" && p.threads == 1)
            strongBase = p.meanTimeSec;
        if (p.study == "weak" && p.threads == 1)
            weakBase = p.meanTimeSec;
    }

    for (auto& p : points) {
        if (p.study == "strong") {
            p.speedup = strongBase / p.meanTimeSec;
            p.efficiency = p.speedup / std::max(1, p.threads);
        } else {
            p.speedup = (weakBase * p.threads) / p.meanTimeSec;
            p.efficiency = weakBase / p.meanTimeSec;
        }
    }

    std::ofstream csv(CSV_FILE.c_str());
    csv << "study,threads,paths,mean_time_sec,stddev_time_sec,speedup,efficiency,price,price_stderr\n";
    for (const auto& p : points) {
        csv << p.study << ","
            << p.threads << ","
            << p.paths << ","
            << p.meanTimeSec << ","
            << p.stddevTimeSec << ","
            << p.speedup << ","
            << p.efficiency << ","
            << p.price << ","
            << p.priceStdErr << "\n";
    }
    csv.close();

    std::cout << "\nSummary\n";
    for (const auto& p : points) {
        std::cout << std::setw(6) << p.study
                  << " threads=" << std::setw(2) << p.threads
                  << " speedup=" << std::setw(8) << p.speedup
                  << " efficiency=" << std::setw(8) << p.efficiency
                  << " time=" << std::setw(8) << p.meanTimeSec << "s\n";
    }
    std::cout << "\nWrote: " << CSV_FILE
              << "\nRun: python " << PLOT_SCRIPT << "\n";
    return 0;
#endif
}
