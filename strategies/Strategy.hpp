#pragma once
#include <vector>

struct BacktestResult {
    double meanPnL    = 0.0;
    double stdPnL     = 0.0;
    double sharpe     = 0.0;
    double maxDrawdown = 0.0;
    double var95      = 0.0;
    double es95       = 0.0;
    std::vector<double> pnlPath; // cumulative PnL series
};

class Strategy {
public:
    virtual ~Strategy() = default;

    virtual BacktestResult run(const std::vector<std::vector<double>>& paths,
                               double rate,
                               double maturity) = 0;

    virtual const char* name() const = 0;
};