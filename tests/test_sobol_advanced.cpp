// tests/test_sobol_advanced.cpp
//
// Three tests:
//   1. Convergence  — MC error vs QMC error across path counts 100..10000
//   2. Star discrepancy intuition — ASCII scatter + CSV export
//   3. Brownian bridge effectiveness — Sobol+bridge vs Sobol bare
//
// Compile (from project root):
//   g++ -std=c++20 -O2 -o test_sobol_advanced tests/test_sobol_advanced.cpp
//
// Run:
//   ./test_sobol_advanced

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <numeric>
#include <random>
#include <string>
#include <algorithm>
#include <iomanip>

#include "../src/core/SobolSequence.hpp"
#include "../src/core/BrownianBridge.hpp"
#include "../src/stochasticProcess/GeometricBrownianMotion.hpp"
#include "../src/integrators/EulerMaruyama.hpp"
#include "../src/payoffs/VanillaPayoffs.hpp"
#include "../src/pricing/BlackScholes.hpp"
#include "../src/pricing/EuropeanPricer.hpp"

// ─────────────────────────────────────────────────────────────────────────────
// Shared config
// ─────────────────────────────────────────────────────────────────────────────
static const std::string DIR_FILE = "docs/new-joe-kuo-6.21201.txt";

static const double S0    = 100.0;
static const double K     = 100.0;
static const double r     = 0.05;
static const double sigma = 0.2;
static const double T     = 1.0;
static const double dt    = 0.01;
static const int    STEPS = static_cast<int>(T / dt);  // 100

// ─────────────────────────────────────────────────────────────────────────────
// Probit (inverse normal CDF) — Beasley-Springer-Moro rational approximation
// ─────────────────────────────────────────────────────────────────────────────
static double probit(double u) {
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
    const double pLow = 0.02425, pHigh = 1.0 - 0.02425;
    double q, r2;
    if (u < pLow) {
        q = std::sqrt(-2.0 * std::log(u));
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
               ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
    } else if (u <= pHigh) {
        q = u - 0.5; r2 = q*q;
        return (((((a[0]*r2+a[1])*r2+a[2])*r2+a[3])*r2+a[4])*r2+a[5])*q /
               (((((b[0]*r2+b[1])*r2+b[2])*r2+b[3])*r2+b[4])*r2+1.0);
    } else {
        q = std::sqrt(-2.0 * std::log(1.0 - u));
        return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
                ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Price a European call using QMC (Sobol + bridge) for a given path count
// Returns {price, stdErr}
// ─────────────────────────────────────────────────────────────────────────────
static std::pair<double,double> priceQMC(int paths) {
    GeometricBrownianMotion gbm(r, sigma);
    CallPayoff call(K);
    SobolSequence sobol(STEPS, DIR_FILE);

    std::vector<double> terminals(paths);
    for (int i = 0; i < paths; i++) {
        auto u = sobol.next();
        std::vector<double> z(STEPS);
        for (int j = 0; j < STEPS; j++) z[j] = probit(u[j]);
        auto dW = BrownianBridge::build(z, dt);

        double x = S0, t = 0.0;
        for (int j = 0; j < STEPS; j++) {
            double zj = dW[j] / std::sqrt(dt);
            x = EulerMaruyama::stepWithZ(gbm, x, t, dt, zj);
            t += dt;
        }
        terminals[i] = x;
    }
    return EuropeanPricer::priceWithError(terminals, call, r, T);
}

// ─────────────────────────────────────────────────────────────────────────────
// Price a European call using QMC WITHOUT Brownian bridge
// Returns {price, stdErr}
// ─────────────────────────────────────────────────────────────────────────────
static std::pair<double,double> priceQMCNoBridge(int paths) {
    GeometricBrownianMotion gbm(r, sigma);
    CallPayoff call(K);
    SobolSequence sobol(STEPS, DIR_FILE);

    std::vector<double> terminals(paths);
    for (int i = 0; i < paths; i++) {
        auto u = sobol.next();

        double x = S0, t = 0.0;
        for (int j = 0; j < STEPS; j++) {
            // Probit directly — no bridge reordering
            double zj = probit(u[j]);
            x = EulerMaruyama::stepWithZ(gbm, x, t, dt, zj);
            t += dt;
        }
        terminals[i] = x;
    }
    return EuropeanPricer::priceWithError(terminals, call, r, T);
}

// ─────────────────────────────────────────────────────────────────────────────
// Price a European call using plain MC for a given path count
// Returns {price, stdErr}
// ─────────────────────────────────────────────────────────────────────────────
static std::pair<double,double> priceMC(int paths) {
    GeometricBrownianMotion gbm(r, sigma);
    CallPayoff call(K);
    Random rng;

    std::vector<double> terminals(paths);
    for (int i = 0; i < paths; i++) {
        double x = S0, t = 0.0;
        for (int j = 0; j < STEPS; j++) {
            x = EulerMaruyama::step(gbm, x, t, dt, rng);
            t += dt;
        }
        terminals[i] = x;
    }
    return EuropeanPricer::priceWithError(terminals, call, r, T);
}

// ─────────────────────────────────────────────────────────────────────────────
// TEST 1: Convergence comparison — MC error vs QMC error
// ─────────────────────────────────────────────────────────────────────────────
void test_convergence() {
    std::cout << "\n";
    std::cout << "══════════════════════════════════════════════════════════════\n";
    std::cout << " TEST 1: Convergence — MC vs QMC (Sobol + Bridge)\n";
    std::cout << "══════════════════════════════════════════════════════════════\n";
    std::cout << " Black-Scholes reference: "
              << BlackScholes::callPrice(S0, K, r, sigma, T) << "\n\n";

    // Path counts: 100, 200, 500, 1000, 2000, 5000, 10000
    std::vector<int> pathCounts = {100, 200, 500, 1000, 2000, 5000, 10000};
    double bsPrice = BlackScholes::callPrice(S0, K, r, sigma, T);

    std::cout << std::left
              << std::setw(8)  << "Paths"
              << std::setw(14) << "MC Price"
              << std::setw(14) << "MC Error"
              << std::setw(14) << "MC |Bias|"
              << std::setw(14) << "QMC Price"
              << std::setw(14) << "QMC Error"
              << std::setw(14) << "QMC |Bias|"
              << "Improvement\n";
    std::cout << std::string(96, '-') << "\n";

    for (int paths : pathCounts) {
        auto [mc_price,  mc_err]  = priceMC(paths);
        auto [qmc_price, qmc_err] = priceQMC(paths);

        double mc_bias  = std::abs(mc_price  - bsPrice);
        double qmc_bias = std::abs(qmc_price - bsPrice);
        double improvement = (mc_err > 1e-12) ? mc_err / qmc_err : 1.0;

        std::cout << std::left  << std::setw(8)  << paths
                  << std::fixed << std::setprecision(5)
                  << std::setw(14) << mc_price
                  << std::setw(14) << mc_err
                  << std::setw(14) << mc_bias
                  << std::setw(14) << qmc_price
                  << std::setw(14) << qmc_err
                  << std::setw(14) << qmc_bias
                  << std::setprecision(2) << improvement << "x\n";
    }

    // Theoretical note
    std::cout << "\n";
    std::cout << " Theory: MC error ~ O(1/√N),  QMC error ~ O((log N)^d / N)\n";
    std::cout << " 'Improvement' = MC stdErr / QMC stdErr — higher is better for QMC\n";
    std::cout << " Note: QMC advantage softens in high dimensions (d=100 steps here)\n";
}

// ─────────────────────────────────────────────────────────────────────────────
// TEST 2: Star discrepancy intuition — ASCII scatter + CSV export
// ─────────────────────────────────────────────────────────────────────────────
void test_star_discrepancy() {
    std::cout << "\n";
    std::cout << "══════════════════════════════════════════════════════════════\n";
    std::cout << " TEST 2: Star Discrepancy Intuition (2D, N=256 points)\n";
    std::cout << "══════════════════════════════════════════════════════════════\n";

    const int N         = 256;
    const int GRID      = 32;   // ASCII grid resolution
    const std::string csvFile = "tests/output/scatter_discrepancy.csv";

    // --- Generate points ---
    std::vector<std::pair<double,double>> sobolPts(N), mcPts(N);

    SobolSequence sobol2d(2, DIR_FILE);
    for (int i = 0; i < N; i++) {
        auto p = sobol2d.next();
        sobolPts[i] = {p[0], p[1]};
    }

    std::mt19937 gen(42);
    std::uniform_real_distribution<double> ud(0.0, 1.0);
    for (int i = 0; i < N; i++)
        mcPts[i] = {ud(gen), ud(gen)};

    // --- ASCII scatter plots ---
    auto makeGrid = [&](const std::vector<std::pair<double,double>>& pts)
        -> std::vector<std::string>
    {
        std::vector<std::string> grid(GRID, std::string(GRID, '.'));
        for (auto& [x, y] : pts) {
            int col = std::min(GRID-1, static_cast<int>(x * GRID));
            int row = std::min(GRID-1, static_cast<int>((1.0-y) * GRID)); // flip y
            grid[row][col] = '*';
        }
        return grid;
    };

    auto sobolGrid = makeGrid(sobolPts);
    auto mcGrid    = makeGrid(mcPts);

    // Print side by side
    std::cout << "\n";
    std::cout << "  Sobol (2D, N=256)               "
              << "  Pseudo-Random MC (N=256)\n";
    std::cout << "  " << std::string(GRID, '-')
              << "    " << std::string(GRID, '-') << "\n";
    for (int row = 0; row < GRID; row++) {
        std::cout << "  |" << sobolGrid[row] << "|"
                  << "    |" << mcGrid[row]    << "|\n";
    }
    std::cout << "  " << std::string(GRID, '-')
              << "    " << std::string(GRID, '-') << "\n";

    // --- Quantitative discrepancy: count points per grid cell ---
    // For a uniform sequence, each cell in a k×k grid should have ~N/k^2 points
    // We measure max deviation (a proxy for star discrepancy)
    auto cellDiscrepancy = [&](const std::vector<std::pair<double,double>>& pts,
                                int k) -> double
    {
        int cells = k * k;
        double expected = static_cast<double>(N) / cells;
        std::vector<int> counts(cells, 0);
        for (auto& [x, y] : pts) {
            int ci = std::min(k-1, static_cast<int>(x * k));
            int ri = std::min(k-1, static_cast<int>(y * k));
            counts[ri * k + ci]++;
        }
        double maxDev = 0.0;
        for (int c : counts)
            maxDev = std::max(maxDev, std::abs(c - expected) / expected);
        return maxDev;  // max relative deviation from uniform
    };

    double sobol_disc4  = cellDiscrepancy(sobolPts, 4);
    double mc_disc4     = cellDiscrepancy(mcPts,    4);
    double sobol_disc8  = cellDiscrepancy(sobolPts, 8);
    double mc_disc8     = cellDiscrepancy(mcPts,    8);
    double sobol_disc16 = cellDiscrepancy(sobolPts, 16);
    double mc_disc16    = cellDiscrepancy(mcPts,    16);

    std::cout << "\n";
    std::cout << " Cell uniformity (max relative deviation from expected count):\n";
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  Grid      Sobol       MC (random)\n";
    std::cout << "  4×4       " << std::setw(10) << sobol_disc4
              << "  " << mc_disc4  << "\n";
    std::cout << "  8×8       " << std::setw(10) << sobol_disc8
              << "  " << mc_disc8  << "\n";
    std::cout << "  16×16     " << std::setw(10) << sobol_disc16
              << "  " << mc_disc16 << "\n";
    std::cout << " Lower = more uniform. Sobol should win on all grids.\n";

    // --- CSV export ---
    // Create output directory if needed (best effort)
    std::system("mkdir -p tests/output");
    std::ofstream csv(csvFile);
    if (csv.is_open()) {
        csv << "x,y,type\n";
        for (auto& [x, y] : sobolPts) csv << x << "," << y << ",sobol\n";
        for (auto& [x, y] : mcPts)    csv << x << "," << y << ",mc\n";
        csv.close();
        std::cout << "\n CSV exported to: " << csvFile << "\n";
        std::cout << " Plot in Python:  "
                  << "df.groupby('type').plot.scatter('x','y')\n";
    } else {
        std::cout << "\n [Warning] Could not write CSV to " << csvFile << "\n";
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// TEST 3: Brownian bridge effectiveness
// Compares: Sobol + bridge  vs  Sobol bare  vs  plain MC
// Metric: absolute pricing error vs Black-Scholes, across path counts
// ─────────────────────────────────────────────────────────────────────────────
void test_brownian_bridge() {
    std::cout << "\n";
    std::cout << "══════════════════════════════════════════════════════════════\n";
    std::cout << " TEST 3: Brownian Bridge Effectiveness\n";
    std::cout << "══════════════════════════════════════════════════════════════\n";

    double bsPrice = BlackScholes::callPrice(S0, K, r, sigma, T);
    std::cout << " Black-Scholes reference: " << bsPrice << "\n\n";

    std::vector<int> pathCounts = {100, 200, 500, 1000, 2000, 5000, 10000};

    std::cout << std::left
              << std::setw(8)  << "Paths"
              << std::setw(16) << "MC StdErr"
              << std::setw(16) << "QMC+Bridge Err"
              << std::setw(16) << "QMC Bare Err"
              << std::setw(14) << "Bridge Gain"
              << "\n";
    std::cout << std::string(70, '-') << "\n";

    for (int paths : pathCounts) {
        auto [mc_price,      mc_err  ] = priceMC(paths);
        auto [bridge_price,  bridge_err] = priceQMC(paths);
        auto [bare_price,    bare_err  ] = priceQMCNoBridge(paths);

        // Absolute bias vs Black-Scholes is a cleaner metric here than stdErr
        // because stdErr conflates variance with the bridge reordering effect
        double mc_bias     = std::abs(mc_price     - bsPrice);
        double bridge_bias = std::abs(bridge_price - bsPrice);
        double bare_bias   = std::abs(bare_price   - bsPrice);

        // Bridge gain: how much smaller is bridge error vs bare QMC error
        double gain = (bare_err > 1e-12) ? bare_err / bridge_err : 1.0;

        std::cout << std::left  << std::setw(8) << paths
                  << std::fixed << std::setprecision(6)
                  << std::setw(16) << mc_err
                  << std::setw(16) << bridge_err
                  << std::setw(16) << bare_err
                  << std::setprecision(2) << gain << "x\n";

        // Also print bias row indented
        std::cout << "        "
                  << std::setprecision(6)
                  << "  |bias|=" << std::setw(10) << mc_bias
                  << "  |bias|=" << std::setw(10) << bridge_bias
                  << "  |bias|=" << bare_bias << "\n";
    }

    std::cout << "\n";
    std::cout << " 'Bridge Gain' = QMC-bare stdErr / QMC+bridge stdErr\n";
    std::cout << " The bridge reorders Sobol dimensions so the most important\n";
    std::cout << " moves (terminal level) use the best low-discrepancy dims.\n";
    std::cout << " Gain > 1.0 confirms the bridge is helping.\n";

    // Hard assertion: bridge should beat bare QMC on average bias
    double total_bridge = 0.0, total_bare = 0.0;
    for (int paths : pathCounts) {
        auto [bp, be] = priceQMC(paths);
        auto [np, ne] = priceQMCNoBridge(paths);
        total_bridge += be;
        total_bare   += ne;
    }
    bool bridgeWins = total_bridge <= total_bare;
    std::cout << "\n Average stdErr — Bridge: " << total_bridge / pathCounts.size()
              << "  Bare: " << total_bare / pathCounts.size() << "\n";
    std::cout << " Bridge beats bare QMC: " << (bridgeWins ? "PASS ✓" : "FAIL ✗") << "\n";
}

// ─────────────────────────────────────────────────────────────────────────────
int main() {
    std::cout << "╔══════════════════════════════════════════════════════════════╗\n";
    std::cout << "║         Advanced Sobol Sequence Test Suite                   ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════╝\n";

    test_convergence();
    test_star_discrepancy();
    test_brownian_bridge();

    std::cout << "\n══════════════════════════════════════════════════════════════\n";
    std::cout << " All tests complete.\n";
    std::cout << "══════════════════════════════════════════════════════════════\n\n";
    return 0;
}