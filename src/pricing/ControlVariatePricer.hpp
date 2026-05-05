// src/pricing/ControlVariatePricer.hpp
#pragma once
#include <vector>
#include <cmath>
#include <numeric>
#include "../payoffs/Payoff.hpp"

class ControlVariatePricer {
public:
    // prices[i]   = terminal asset price for path i
    // payoffs[i]  = payoff(prices[i]) — pre-computed
    // Analytical mean of the control (E[S_T] = S0 * exp(r*T))
    static double price(
        const std::vector<double>& terminalPrices,
        const Payoff& payoff,
        double r,
        double T,
        double S0)
    {
        int n = terminalPrices.size();
        double disc        = std::exp(-r * T);
        double analyticalMean = S0 * std::exp(r * T);  // E[S_T] under risk-neutral

        // Compute raw payoffs and control values
        std::vector<double> Y(n), Z(n);
        for (int i = 0; i < n; i++) {
            Y[i] = payoff(terminalPrices[i]);
            Z[i] = terminalPrices[i];
        }

        double meanY = std::accumulate(Y.begin(), Y.end(), 0.0) / n;
        double meanZ = std::accumulate(Z.begin(), Z.end(), 0.0) / n;

        // Estimate optimal beta: cov(Y,Z) / var(Z)
        double covYZ = 0.0, varZ = 0.0;
        for (int i = 0; i < n; i++) {
            covYZ += (Y[i] - meanY) * (Z[i] - meanZ);
            varZ  += (Z[i] - meanZ) * (Z[i] - meanZ);
        }
        double beta = (varZ > 1e-12) ? covYZ / varZ : 0.0;

        // Adjusted estimator: Y - beta*(Z - E[Z])
        double sum = 0.0;
        for (int i = 0; i < n; i++)
            sum += Y[i] - beta * (Z[i] - analyticalMean);

        return disc * sum / n;
    }

    // Returns {price, variance_reduction_ratio} for diagnostics
    static std::pair<double, double> priceWithDiagnostics(
        const std::vector<double>& terminalPrices,
        const Payoff& payoff,
        double r,
        double T,
        double S0)
    {
        int n = terminalPrices.size();
        double disc           = std::exp(-r * T);
        double analyticalMean = S0 * std::exp(r * T);

        std::vector<double> Y(n), Z(n);
        for (int i = 0; i < n; i++) {
            Y[i] = payoff(terminalPrices[i]);
            Z[i] = terminalPrices[i];
        }

        double meanY = std::accumulate(Y.begin(), Y.end(), 0.0) / n;
        double meanZ = std::accumulate(Z.begin(), Z.end(), 0.0) / n;

        double covYZ = 0.0, varZ = 0.0, varY = 0.0;
        for (int i = 0; i < n; i++) {
            covYZ += (Y[i] - meanY) * (Z[i] - meanZ);
            varZ  += (Z[i] - meanZ) * (Z[i] - meanZ);
            varY  += (Y[i] - meanY) * (Y[i] - meanY);
        }
        double beta = (varZ > 1e-12) ? covYZ / varZ : 0.0;

        // Variance of adjusted estimator: var(Y) - cov^2/var(Z)
        double varAdjusted = varY - covYZ * covYZ / varZ;
        double reductionRatio = (varY > 1e-12) ? varAdjusted / varY : 1.0;

        double sum = 0.0;
        for (int i = 0; i < n; i++)
            sum += Y[i] - beta * (Z[i] - analyticalMean);

        return {disc * sum / n, reductionRatio};
    }
};