// src/simulators/ParallelMonteCarloSimulator.hpp
#pragma once
#include <vector>
#include <omp.h>
#include "../stochasticProcess/StochasticProcess.hpp"
#include "../core/Random.hpp"

template<typename Integrator>
class ParallelMonteCarloSimulator {
public:
    static std::vector<double> simulate(
        const StochasticProcess& process,
        double x0, double T, double dt, int paths,
        int numThreads = 0)  // 0 = let OpenMP decide
    {
        int steps = static_cast<int>(T / dt);
        std::vector<double> results(paths);

        if (numThreads > 0)
            omp_set_num_threads(numThreads);

        // Each thread gets its own RNG seeded differently
        int maxThreads = omp_get_max_threads();
        std::vector<Random> rngs;
        rngs.reserve(maxThreads);
        for (int i = 0; i < maxThreads; i++)
            rngs.emplace_back(42 + i * 1000);  // deterministic but distinct seeds

        #pragma omp parallel for schedule(dynamic, 64)
        for (int i = 0; i < paths; i++) {
            int tid = omp_get_thread_num();
            Random& rng = rngs[tid];

            double x = x0, t = 0.0;
            for (int j = 0; j < steps; j++) {
                x = Integrator::step(process, x, t, dt, rng);
                t += dt;
            }
            results[i] = x;
        }

        return results;
    }
};