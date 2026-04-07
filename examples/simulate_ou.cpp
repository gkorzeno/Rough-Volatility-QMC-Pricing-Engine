#include <iostream>
#include <numeric>   // for std::accumulate
#include <cmath>     // for pow
#include "../src/integrators/EulerMaruyama.hpp"
#include "../src/stochasticProcess/OrnsteinUhlenbeck.hpp"
#include "../src/simulators/MonteCarloSimulator.hpp"

int main() {

    OrnsteinUhlenbeck ou(1.0, 0.0, 0.3);

    double x0 = 2.0;
    double T = 2.0;
    double dt = 0.01;
    int paths = 10000;

    auto results = MonteCarloSimulator<EulerMaruyama>::simulate(ou, x0, T, dt, paths);

    // Compute sample mean
    double sum = std::accumulate(results.begin(), results.end(), 0.0);
    double mean = sum / paths;

    // Compute sample variance
    double sq_sum = 0.0;
    for (double x : results) {
        sq_sum += (x - mean) * (x - mean);
    }
    double variance = sq_sum / paths; // Monte Carlo variance estimate

    // Output
    std::cout << "Simulating: " << paths << " paths" << std::endl;
    std::cout << "Estimated E[X_T] = " << mean << std::endl; 
    std::cout << "Estimated Variance = " << variance << std::endl;

    // Optional: print theoretical mean & variance
    std::cout << "Theoretical E[X_T] = " << ou.theoreticalMean(x0, T) << std::endl;
    std::cout << "Theoretical Variance = " << ou.theoreticalVariance(T) << std::endl;

    return 0;
}
