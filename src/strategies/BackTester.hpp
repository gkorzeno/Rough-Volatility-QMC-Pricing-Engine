#pragma once
#include "../strategies/Strategy.hpp"
#include <iostream>
#include <vector>

class Backtester {
public:
    Backtester(double rate, double maturity)
        : r(rate), T(maturity) {}

    BacktestResult run(Strategy& strategy,
                       const std::vector<std::vector<double>>& paths) const {
        return strategy.run(paths, r, T);
    }

    void report(const Strategy& strategy,
                const BacktestResult& result) const {
        std::cout << "=== Strategy: " << strategy.name() << " ===\n"
                  << "Mean PnL:  " << result.meanPnL    << "\n"
                  << "Std Dev:   " << result.stdPnL     << "\n"
                  << "Sharpe:    " << result.sharpe     << "\n"
                  << "Max DD:    " << result.maxDrawdown << "\n"
                  << "VaR 95%:   " << result.var95      << "\n"
                  << "ES 95%:    " << result.es95       << "\n";
    }

private:
    double r;
    double T;
};