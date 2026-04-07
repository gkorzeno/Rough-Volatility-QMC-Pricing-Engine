#include <iostream>
#include <numeric>
#include "../src/integrators/EulerMaruyama.hpp"
#include "../src/stochasticProcess/GeometricBrownianMotion.hpp"
#include "../src/simulators/MonteCarloSimulator.hpp"

int main() {

    GeometricBrownianMotion gbm(0.05, 0.2);

    double x0 = 100.0;
    double T = 1.0;
    double dt = 0.01;
    int paths = 10000;

    auto results = MonteCarloSimulator<EulerMaruyama>::simulate(gbm, x0, T, dt, paths);

    double mean = std::accumulate(results.begin(), results.end(), 0.0) / paths;

    std::cout << "Estimated E[X_T] = " << mean << std::endl;

}
