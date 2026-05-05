#pragma once
#include "Strategy.hpp"
#include "../pricing/BlackScholes.hpp"
#include <cmath>
#include <algorithm>

class ShortStraddleStrategy : public Strategy {
public:
    ShortStraddleStrategy(double strike,
                          double impliedVol,
                          double transactionCost = 0.0)
        : K(strike), sigma(impliedVol), tc(transactionCost) {}

    BacktestResult run(const std::vector<std::vector<double>>& paths,
                       double r,
                       double T) override;

    const char* name() const override { return "Short Straddle"; }

private:
    double K;     // strike
    double sigma; // implied vol used to price the straddle at inception
    double tc;    // one-way transaction cost (subtracted from premium)
};