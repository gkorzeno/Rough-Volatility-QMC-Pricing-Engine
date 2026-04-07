// examples/simulate_heston.cpp
#include <iostream>
#include <numeric>
#include <cmath>
#include "../src/stochasticProcess/Heston.hpp"
#include "../src/simulators/MultiDimensionalMonteCarloSimulator.hpp"

int main() {
    // Standard Heston parameters
    double mu    = 0.05;
    double kappa = 2.0;
    double vbar  = 0.04;   // long-run variance = 20% vol
    double xi    = 0.3;    // vol of vol
    double rho   = -0.7;   // typical negative correlation

    Heston heston(mu, kappa, vbar, xi, rho);

    std::cout << "Feller condition satisfied: "
              << (heston.fellerConditionSatisfied() ? "yes" : "no") << std::endl;

    std::vector<double> x0 = {100.0, 0.04};  // S0=100, v0=4%
    double T  = 1.0;
    double dt = 0.005;    // finer dt — Heston needs it for variance stability
    int paths = 10000;

    std::cout << "diffusionIsCholesky: " << heston.diffusionIsCholesky() << std::endl;

    // Also manually test one step
    Random rng;
    auto test_drift = heston.drift(x0, 0.0);
    auto test_diff  = heston.diffusion(x0, 0.0);
    std::cout << "Drift: " << test_drift[0] << ", " << test_drift[1] << std::endl;
    std::cout << "Diffusion L[0][0]: " << test_diff[0][0] << std::endl;

    auto results = MultiDimensionalMonteCarloSimulator::simulate(
        heston, x0, T, dt, paths);

    // Collect terminal S and v
    double sumS = 0.0, sumV = 0.0;
    int neg_var = 0;
    for (auto& r : results) {
        sumS += r[0];
        sumV += r[1];
        if (r[1] < 0) neg_var++;
    }

    double meanS = sumS / paths;
    double meanV = sumV / paths;

    std::cout << "E[S_T] simulated  = " << meanS << std::endl;
    std::cout << "E[S_T] theoretical = " << x0[0] * std::exp(mu * T) << std::endl;
    std::cout << "E[v_T] simulated  = " << meanV << std::endl;
    std::cout << "E[v_T] theoretical = " << vbar + (x0[1] - vbar)
                                            * std::exp(-kappa * T) << std::endl;
    std::cout << "Paths with negative variance: " << neg_var << std::endl;

    return 0;
}