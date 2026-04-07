#include <iostream>
#include "../src/stochasticProcess/GeometricBrownianMotion.hpp"
#include "../src/simulators/MonteCarloSimulator.hpp"
#include "../src/integrators/EulerMaruyama.hpp"

int main() {

    GeometricBrownianMotion gbm(0.05, 0.2);

    double S0 = 100;
    double T  = 1.0;
    double dt = 0.01;
    int paths = 100000;

    auto results =
        MonteCarloSimulator<EulerMaruyama>::simulate(
            gbm, S0, T, dt, paths);

    double sum = 0;
    for (double x : results) sum += x;

    double mean = sum / paths;

    std::cout << "Estimated E[S_T] = " << mean << std::endl;
}