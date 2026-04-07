// src/pricing/MultiAssetPricer.hpp
#pragma once
#include <vector>
#include <cmath>
#include <numeric>
#include "../payoffs/MultiAssetPayoff.hpp"

class MultiAssetPricer {
public:
    static std::pair<double, double> price(
        const std::vector<std::vector<double>>& terminalPrices,
        const MultiAssetPayoff& payoff,
        double r,
        double T)
    {
        int n = terminalPrices.size();
        std::vector<double> payoffs(n);
        for (int i = 0; i < n; i++)
            payoffs[i] = payoff(terminalPrices[i]);

        double mean = std::accumulate(payoffs.begin(), payoffs.end(), 0.0) / n;

        double sq_sum = 0.0;
        for (double p : payoffs) sq_sum += (p - mean) * (p - mean);
        double stdErr = std::sqrt(sq_sum / n) / std::sqrt(n);

        double disc = std::exp(-r * T);
        return {disc * mean, disc * stdErr};
    }
};