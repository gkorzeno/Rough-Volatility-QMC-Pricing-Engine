// examples/benchmark_parallel.cpp
#include <iostream>
#include <numeric>
#include "../src/stochasticProcess/GeometricBrownianMotion.hpp"
#include "../src/simulators/MonteCarloSimulator.hpp"
#include "../src/simulators/ParallelMonteCarloSimulator.hpp"
#include "../src/simulators/ParallelAntitheticSimulator.hpp"
#include "../src/integrators/EulerMaruyama.hpp"
#include "../src/payoffs/VanillaPayoffs.hpp"
#include "../src/pricing/EuropeanPricer.hpp"
#include "../src/pricing/BlackScholes.hpp"
#include "../src/utils/Benchmark.hpp"
#include <omp.h>

int main() {
    GeometricBrownianMotion gbm(0.05, 0.2);
    CallPayoff call(100.0);

    double S0 = 100.0, r = 0.05, T = 1.0, dt = 0.01;
    int paths = 200000;

    std::cout << "Threads available: " << omp_get_max_threads() << "\n";
    std::cout << "Black-Scholes: "
              << BlackScholes::callPrice(S0, 100.0, r, 0.2, T) << "\n\n";

    std::vector<double> res_serial, res_parallel, res_parallel_anti;

    Benchmark::run("Serial   (200k paths)", [&]() {
        res_serial = MonteCarloSimulator<EulerMaruyama>::simulate(
            gbm, S0, T, dt, paths);
    });

    Benchmark::run("Parallel (200k paths)", [&]() {
        res_parallel = ParallelMonteCarloSimulator<EulerMaruyama>::simulate(
            gbm, S0, T, dt, paths);
    });

    Benchmark::run("Parallel antithetic (200k paths)", [&]() {
        res_parallel_anti = ParallelAntitheticSimulator<EulerMaruyama>::simulate(
            gbm, S0, T, dt, paths);
    });

    auto [p1, e1] = EuropeanPricer::priceWithError(res_serial,        call, r, T);
    auto [p2, e2] = EuropeanPricer::priceWithError(res_parallel,      call, r, T);
    auto [p3, e3] = EuropeanPricer::priceWithError(res_parallel_anti, call, r, T);

    std::cout << "\nSerial price:             " << p1 << " ± " << e1 << "\n";
    std::cout << "Parallel price:           " << p2 << " ± " << e2 << "\n";
    std::cout << "Parallel antithetic:      " << p3 << " ± " << e3 << "\n";

    // Thread scaling test
    std::cout << "\n-- Thread scaling (200k paths, Euler-Maruyama) --\n";
    for (int t : {1, 2, 4, 8}) {
        if (t > omp_get_max_threads()) break;
        Benchmark::run("  " + std::to_string(t) + " thread(s)", [&]() {
            ParallelMonteCarloSimulator<EulerMaruyama>::simulate(
                gbm, S0, T, dt, paths, t);
        });
    }

    return 0;
}