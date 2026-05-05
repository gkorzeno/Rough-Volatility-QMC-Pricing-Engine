// src/pricing/HestonPricer.hpp
#pragma once
#include <vector>
#include <cmath>
#include <numeric>
#include "../payoffs/Payoff.hpp"
#include "../stochasticProcess/Heston.hpp"
#include "../simulators/MultiDimensionalMonteCarloSimulator.hpp"

class HestonPricer {
public:
    static std::pair<double, double> priceEuropean(
        const Heston& heston,
        const Payoff& payoff,
        double S0, double v0,
        double r, double T, double dt,
        int paths)
    {
        std::vector<double> x0 = {S0, v0};
        auto results = MultiDimensionalMonteCarloSimulator::simulate(
            heston, x0, T, dt, paths);

        int n = results.size();
        std::vector<double> payoffs(n);
        for (int i = 0; i < n; i++)
            payoffs[i] = payoff(results[i][0]);  // x[0] = S

        double mean = std::accumulate(payoffs.begin(), payoffs.end(), 0.0) / n;

        double sq_sum = 0.0;
        for (double p : payoffs) sq_sum += (p - mean) * (p - mean);
        double stdErr = std::sqrt(sq_sum / n) / std::sqrt(n);

        double disc = std::exp(-r * T);
        return {disc * mean, disc * stdErr};
    }
};