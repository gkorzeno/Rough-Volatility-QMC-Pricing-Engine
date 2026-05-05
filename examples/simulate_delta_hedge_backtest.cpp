#include <iostream>
#include <iomanip>
#include "../src/stochasticProcess/GeometricBrownianMotion.hpp"
#include "../src/simulators/PathSimulator.hpp"
#include "../src/integrators/EulerMaruyama.hpp"
#include "../src/strategies/ShortStraddleStrategy.hpp"
#include "../src/strategies/DeltaHedgedStraddleStrategy.hpp"
#include "../src/strategies/Backtester.hpp"

static void printRow(const char* label, const BacktestResult& r) {
    std::cout << std::left  << std::setw(34) << label
              << std::right << std::fixed << std::setprecision(4)
              << std::setw(10) << r.meanPnL
              << std::setw(10) << r.stdPnL
              << std::setw(10) << r.sharpe
              << std::setw(10) << r.maxDrawdown
              << std::setw(10) << r.var95
              << std::setw(10) << r.es95 << "\n";
}

int main() {
    const double S0          = 100.0;
    const double K           = 100.0;
    const double r           = 0.05;
    const double realisedVol = 0.20;   // true GBM vol
    const double impliedVol  = 0.25;   // we sell the straddle at a 5-vol premium
    const double T           = 1.0;
    const double dt          = 1.0 / 252; // daily path steps
    const int    paths       = 20000;
    const double tc          = 0.001;

    GeometricBrownianMotion gbm(r, realisedVol);
    auto allPaths = PathSimulator<EulerMaruyama>::simulate(gbm, S0, T, dt, paths);

    Backtester bt(r, T);

    // Unhedged: keep the raw short straddle to expiry
    ShortStraddleStrategy raw(K, impliedVol, tc);

    // Hedged variants at different rebalancing frequencies
    DeltaHedgedStraddleStrategy daily(K, impliedVol, 1.0/252, tc);   // daily
    DeltaHedgedStraddleStrategy weekly(K, impliedVol, 1.0/52,  tc);  // weekly
    DeltaHedgedStraddleStrategy monthly(K, impliedVol, 1.0/12, tc);  // monthly

    const auto rawResult     = bt.run(raw,     allPaths);
    const auto dailyResult   = bt.run(daily,   allPaths);
    const auto weeklyResult  = bt.run(weekly,  allPaths);
    const auto monthlyResult = bt.run(monthly, allPaths);

    std::cout << "\nRealised vol: " << realisedVol
              << "  |  Implied vol sold at: " << impliedVol
              << "  |  Paths: " << paths << "\n\n";

    std::cout << std::left  << std::setw(34) << "Strategy"
              << std::right
              << std::setw(10) << "MeanPnL"
              << std::setw(10) << "StdPnL"
              << std::setw(10) << "Sharpe"
              << std::setw(10) << "MaxDD"
              << std::setw(10) << "VaR95"
              << std::setw(10) << "ES95" << "\n";
    std::cout << std::string(84, '-') << "\n";

    printRow("Short Straddle (unhedged)",         rawResult);
    printRow("Delta-Hedged (daily rebalance)",     dailyResult);
    printRow("Delta-Hedged (weekly rebalance)",    weeklyResult);
    printRow("Delta-Hedged (monthly rebalance)",   monthlyResult);

    std::cout << "\n";
    bt.report(daily, dailyResult);
    return 0;
}