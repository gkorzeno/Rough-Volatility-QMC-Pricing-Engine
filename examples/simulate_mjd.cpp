// tests/testMertonJumpDiffusion.cpp
#include <iostream>
#include <numeric>
#include <cmath>
#include <vector>
#include "../src/stochasticProcess/MertonJumpDiffusion.hpp"
#include "../src/simulators/MonteCarloSimulator.hpp"

int main() {
    // Parameters: drift, diffusion, jump rate, mean jump size, jump vol
    // lambda=1.0 means ~1 jump per year on average — visible but not overwhelming
    double mu     = 0.05;
    double sigma  = 0.2;
    double lambda = 1.0;
    double muJ    = -0.1;   // negative mean: jumps tend downward (crash model)
    double sigmaJ = 0.15;

    MertonJumpDiffusion mjd(mu, sigma, lambda, muJ, sigmaJ);

    double x0    = 100.0;
    double T     = 1.0;
    double dt    = 0.01;
    int    paths = 50000;   // more paths than GBM test — jumps add variance

    auto results = MonteCarloSimulator<EulerMaruyama>::simulate(mjd, x0, T, dt, paths);

    // --- Basic statistics ---
    double sum = std::accumulate(results.begin(), results.end(), 0.0);
    double mean = sum / paths;

    double sq_sum = 0.0;
    for (double x : results) sq_sum += (x - mean) * (x - mean);
    double variance = sq_sum / paths;

    // --- Jump verification: count paths that moved "too far" for pure diffusion ---
    // Pure GBM 3-sigma move in 1yr ~ x0 * exp((mu - 0.5*sigma^2)*T ± 3*sigma*sqrt(T))
    double diffusion_only_3sigma = x0 * std::exp((mu - 0.5*sigma*sigma)*T 
                                                  + 3*sigma*std::sqrt(T));
    int extreme_moves = 0;
    double min_val = results[0], max_val = results[0];
    for (double x : results) {
        if (x > diffusion_only_3sigma) extreme_moves++;
        if (x < min_val) min_val = x;
        if (x > max_val) max_val = x;
    }

    // --- Output ---
    std::cout << "=== Merton Jump Diffusion Test ===" << std::endl;
    std::cout << "Paths: " << paths << ", T=" << T << ", dt=" << dt << std::endl;
    std::cout << std::endl;

    std::cout << "-- Mean --" << std::endl;
    std::cout << "Simulated  E[X_T] = " << mean << std::endl;
    std::cout << "Theoretical E[X_T] = " << mjd.theoreticalMean(x0, T) << std::endl;
    std::cout << std::endl;

    std::cout << "-- Variance --" << std::endl;
    std::cout << "Simulated  Var[X_T] = " << variance << std::endl;
    std::cout << "Theoretical Var[X_T] = " << mjd.theoreticalVariance(x0, T) << std::endl;
    std::cout << std::endl;

    std::cout << "-- Jump Diagnostics --" << std::endl;
    std::cout << "Min path value: " << min_val << std::endl;
    std::cout << "Max path value: " << max_val << std::endl;
    std::cout << "Paths beyond pure-diffusion 3-sigma: " << extreme_moves 
              << " (" << 100.0 * extreme_moves / paths << "%)" << std::endl;
    std::cout << "(Expected ~few % for lambda=1 with positive jump component)" << std::endl;

    return 0;
}