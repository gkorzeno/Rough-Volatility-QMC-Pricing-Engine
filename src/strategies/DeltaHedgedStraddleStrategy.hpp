#pragma once
#include "Strategy.hpp"
#include "../pricing/BlackScholes.hpp"
#include <cmath>
#include <stdexcept>

class DeltaHedgedStraddleStrategy : public Strategy {
public:
    DeltaHedgedStraddleStrategy(double strike,
                                double impliedVol,
                                double hedgeInterval,
                                double transactionCost = 0.0)
        : K(strike), sigma(impliedVol), dtHedge(hedgeInterval), tc(transactionCost)
    {
        if (dtHedge <= 0.0)
            throw std::invalid_argument("hedgeInterval must be positive");
    }

    BacktestResult run(const std::vector<std::vector<double>>& paths,
                       double r, double T) override;

    const char* name() const override { return "Delta-Hedged Short Straddle"; }

private:
    double K, sigma, dtHedge, tc;

    // Delta of a LONG straddle = delta_call + delta_put = N(d1) + (N(d1)-1) = 2N(d1)-1.
    // We are SHORT, so our delta position is -(2N(d1)-1).
    // Near ATM this is ~0, which is correct: a straddle has near-zero delta at the money.
    static double straddleDelta(double S, double K, double r,
                                double sigma, double tau) {
        if (tau <= 0.0) return 0.0;
        const double sqrtTau = std::sqrt(tau);
        const double d1 = (std::log(S / K) + (r + 0.5 * sigma * sigma) * tau)
                        / (sigma * sqrtTau);
        const double Nd1 = 0.5 * std::erfc(-d1 / std::sqrt(2.0));
        return 2.0 * Nd1 - 1.0;  // long straddle delta; caller negates for short
    }

    static double straddlePrice(double S, double K, double r,
                                double sigma, double tau) {
        return BlackScholes::callPrice(S, K, r, sigma, tau)
             + BlackScholes::putPrice(S, K, r, sigma, tau);
    }
};