#include <iostream>
#include <iomanip>
#include "../src/stochasticProcess/GeometricBrownianMotion.hpp"
#include "../src/simulators/PathSimulator.hpp"
#include "../src/integrators/EulerMaruyama.hpp"
#include "../src/strategies/Backtester.hpp"
#include "../src/strategies/vrp/VRPStrategy.hpp"
#include "../src/strategies/mispricing/MispricingStrategy.hpp"
#include "../src/strategies/dispersion/DispersionStrategy.hpp"
#include "../src/strategies/termstructure/TermStructureStrategy.hpp"
#include "../src/roughvol/RoughBergomiModel.hpp"
#include "../src/stochasticProcess/CorrelatedGBM.hpp"

static void printRow(const char* label, const BacktestResult& r) {
    std::cout << std::left  << std::setw(36) << label
              << std::right << std::fixed << std::setprecision(4)
              << std::setw(10) << r.meanPnL
              << std::setw(10) << r.stdPnL
              << std::setw(10) << r.sharpe
              << std::setw(10) << r.maxDrawdown
              << std::setw(10) << r.var95 << "\n";
}

int main() {
    const double S0     = 100.0;
    const double r      = 0.05;
    const double T      = 1.0;
    const double dt     = 1.0 / 252;
    const double sigmaR = 0.20;
    const double sigmaI = 0.25;
    const double K      = 100.0;
    const int    nPaths = 10000;
    const double tc     = 0.001;

    GeometricBrownianMotion gbm(r, sigmaR);
    auto paths = PathSimulator<EulerMaruyama>::simulate(gbm, S0, T, dt, nPaths);

    RoughBergomiModel::Parameters rbParams{0.04, 1.5, -0.7, 0.12};

    // ── Build a realistic surface with proper skew and term structure ──────
    // Use a log-moneyness parameterisation so skew is meaningful.
    // Market: steep put skew (equity-style), upward sloping term structure.
    const std::vector<double> strikes  = {80, 85, 90, 95, 100, 105, 110, 115, 120};
    const std::vector<double> expiries = {0.25, 0.5, 1.0, 2.0};

    // Market surface: strong put skew + vol term structure slope of ~0.03/year
    std::vector<std::vector<double>> marketVols(expiries.size(),
        std::vector<double>(strikes.size()));
    for (std::size_t ti = 0; ti < expiries.size(); ++ti) {
        for (std::size_t ki = 0; ki < strikes.size(); ++ki) {
            const double k    = std::log(strikes[ki] / S0);
            const double base = 0.20
                              - 0.12 * k          // strong put skew
                              + 0.08 * k * k;     // convexity / smile
            const double termAdj = 0.03 * expiries[ti]; // 3 vol pts/year slope
            marketVols[ti][ki] = std::max(0.05, base + termAdj);
        }
    }
    VolSurface marketSurface(strikes, expiries, marketVols);

    // SVI surface: models the market well but overshoots wings by ~1.5 vol pts
    // → wings look cheap in market, models say they're worth more
    // std::vector<std::vector<double>> sviVols = marketVols;
    // for (std::size_t ti = 0; ti < expiries.size(); ++ti)
    //     for (std::size_t ki = 0; ki < strikes.size(); ++ki) {
    //         const double k = std::log(strikes[ki] / S0);
    //         // Wings (|k| > 0.1) are underpriced in market vs SVI
    //         sviVols[ti][ki] += 0.015 * k * k;
    //     }
    std::vector<std::vector<double>> sviVols = marketVols;
    for (std::size_t ti = 0; ti < expiries.size(); ++ti)
        for (std::size_t ki = 0; ki < strikes.size(); ++ki) {
            const double k = std::log(strikes[ki] / S0);
            // sviVols[ti][ki] += 0.025 * k * k + 0.005 * std::abs(k);
            sviVols[ti][ki] += 0.060 * k * k + 0.015 * std::abs(k);
        }
    VolSurface sviSurface(strikes, expiries, sviVols);

    // rBergomi surface: put skew steeper than market
    // std::vector<std::vector<double>> rbVols = marketVols;
    // for (std::size_t ti = 0; ti < expiries.size(); ++ti)
    //     for (std::size_t ki = 0; ki < strikes.size(); ++ki) {
    //         const double k = std::log(strikes[ki] / S0);
    //         rbVols[ti][ki] += 0.012 * k * k + 0.008 * k; // more skew + convexity
    //     }
    std::vector<std::vector<double>> rbVols = marketVols;
    for (std::size_t ti = 0; ti < expiries.size(); ++ti)
        for (std::size_t ki = 0; ki < strikes.size(); ++ki) {
            const double k = std::log(strikes[ki] / S0);
            // rbVols[ti][ki] += 0.020 * k * k - 0.010 * k; // more put skew
            rbVols[ti][ki] += 0.035 * k * k - 0.030 * k;
        }
    VolSurface rbSurface(strikes, expiries, rbVols);

    // Local vol surface: flatter skew than market (Dupire tends to flatten)
    // std::vector<std::vector<double>> lvVols = marketVols;
    // for (std::size_t ti = 0; ti < expiries.size(); ++ti)
    //     for (std::size_t ki = 0; ki < strikes.size(); ++ki) {
    //         const double k = std::log(strikes[ki] / S0);
    //         lvVols[ti][ki] -= 0.010 * k * k; // less convexity → OTM calls cheap
    //     }
    std::vector<std::vector<double>> lvVols = marketVols;
    for (std::size_t ti = 0; ti < expiries.size(); ++ti)
        for (std::size_t ki = 0; ki < strikes.size(); ++ki) {
            const double k = std::log(strikes[ki] / S0);
            // lvVols[ti][ki] -= 0.018 * k * k; // significantly flatter wings
            lvVols[ti][ki] -= 0.040 * k * k;
        }
    VolSurface lvSurface(strikes, expiries, lvVols);

    Backtester bt(r, T);

    std::cout << "\n";
    std::cout << std::left << std::setw(36) << "Strategy"
              << std::right
              << std::setw(10) << "MeanPnL"
              << std::setw(10) << "StdPnL"
              << std::setw(10) << "Sharpe"
              << std::setw(10) << "MaxDD"
              << std::setw(10) << "VaR95" << "\n"
              << std::string(86, '-') << "\n";

    // ── 1. VRP ────────────────────────────────────────────────────────────
    {
        VRPStrategy::Config cfg{rbParams, sigmaI, K,
            1.0/52, 0.02, tc, 2000, 20};
        VRPStrategy vrp(cfg);
        const auto sig = vrp.computeSignal(S0, r, T);
        std::cout << "VRP  IV=" << std::fixed << std::setprecision(4) << sig.impliedVol
                  << "  forecast=" << sig.forecastedVol
                  << "  VRP=" << sig.vrp
                  << (sig.shouldTrade ? "  TRADE\n" : "  FLAT\n");
        printRow("VRP Strategy", bt.run(vrp, paths));
    }

    // ── 2. Mispricing ─────────────────────────────────────────────────────
    {
        // minMispricing = 0.01 (1 vol point) — achievable with the surfaces above
        MispricingStrategy::Config cfg{
            marketSurface, sviSurface, rbSurface, lvSurface,
            S0, 1.0/52, tc,
            /*minMispricing=*/ 0.003,
            /*maxLegs=*/       3};
        MispricingStrategy misp(cfg);
        auto res = misp.runDetailed(paths, r, T);
        std::cout << "\nMispricing: " << res.mispricings.size()
                  << " opportunities, " << res.trades.size() << " trades\n";
        for (const auto& t : res.trades)
            std::cout << "  " << t.description
                      << "  netEdge=" << std::setprecision(4) << t.netEdge << "\n";
        printRow("Mispricing RV", res.aggregate);

        {
            auto allPoints = MispricingScanner::scan(
            marketSurface, sviSurface, rbSurface, lvSurface, 0.0); // threshold=0
            std::cout << "\nMispricing scan (all points, no threshold):\n";
            for (std::size_t i = 0; i < std::min<std::size_t>(5, allPoints.size()); ++i) {
                const auto& p = allPoints[i];
                std::cout << "  T=" << p.expiry << " K=" << p.strike
                    << " mkt=" << p.marketIV
                    << " svi=" << p.sviIV
                    << " rb=" << p.rbergomiIV
                    << " consensus=" << p.consensusMispricing << "\n";
            }
        }
    }

    // ── 3. Dispersion ─────────────────────────────────────────────────────
    {
        // Realistic dispersion: index IV < avg constituent IV
        // Index vol = 0.22, constituents at 0.25, 0.27, 0.30 → spread = 0.063
        const std::vector<double> mu    = {r, r, r};
        const std::vector<double> sigma = {0.20, 0.22, 0.24};
        const std::vector<std::vector<double>> corr = {
            {1.0, 0.55, 0.40},
            {0.55, 1.0, 0.45},
            {0.40, 0.45, 1.0}};
        CorrelatedGBM cgbm(mu, sigma, corr);

        DispersionStrategy::Config cfg{
            /*indexStrike*/     K,
            /*indexImpliedVol*/ 0.22,
            /*basketWeights*/   {1.0/3, 1.0/3, 1.0/3},
            /*constStrikes*/    {K, K, K},
            /*constIVs*/        {0.25, 0.27, 0.30},
            cgbm,
            /*spotVec*/         {S0, S0, S0},
            /*hedgeInterval*/   1.0/52,
            /*tc*/              tc,
            /*entryThreshold*/  0.03   // 3 vol pt spread → clearly achievable
        };
        DispersionStrategy disp(cfg);
        printRow("Dispersion Strategy", bt.run(disp, paths));
    }

    // ── 4. Term structure / skew ──────────────────────────────────────────
    {
        TermStructureStrategy::Config cfg{
            marketSurface, S0, r,
            /*hedgeInterval*/      1.0/52,
            /*tc*/                 tc,
            /*skewThreshold*/      0.03,   // normalised skew — surface gives ~0.09
            /*termSpreadThresh*/   0.01,   // 1 vol pt spread — surface gives ~0.03
            /*maxSignals*/         3};
        TermStructureStrategy ts(cfg);
        auto res = ts.runDetailed(paths, r, T);

        std::cout << "\nTerm structure signals: " << res.signals.size() << "\n";
        for (const auto& sig : res.signals)
            std::cout << "  ["
                      << (sig.type == SkewAnalyzer::Signal::CalendarSpread
                              ? "Calendar" : "RiskRev")
                      << "] " << sig.description
                      << "  edge=" << std::setprecision(4) << sig.expectedEdge << "\n";

        std::cout << "Skew by expiry:\n";
        for (const auto& sk : res.skewByExpiry)
            std::cout << "  T=" << sk.expiry
                      << " ATM=" << sk.atmVol
                      << " normSkew=" << sk.skew
                      << " RR=" << sk.riskReversalIV << "\n";

        printRow("Term Structure / Skew", res.aggregate);
    }

    std::cout << "\n";
    return 0;
}