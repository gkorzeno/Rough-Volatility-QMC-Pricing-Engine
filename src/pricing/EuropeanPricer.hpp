// src/pricing/EuropeanPricer.hpp
#pragma once
#include <vector>
#include <numeric>
#include <cmath>
#include "../payoffs/Payoff.hpp"

class EuropeanPricer {
public:
    static double price(
        const std::vector<double>& terminalPrices,
        const Payoff& payoff,
        double r,
        double T)
    {
        double sum = 0.0;
        for (double S : terminalPrices)
            sum += payoff(S);

        double mean = sum / terminalPrices.size();
        return std::exp(-r * T) * mean;
    }

    // Also return standard error so you know confidence interval
    static std::pair<double, double> priceWithError(
        const std::vector<double>& terminalPrices,
        const Payoff& payoff,
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

        double discount = std::exp(-r * T);
        return {discount * mean, discount * stdErr};
    }
};