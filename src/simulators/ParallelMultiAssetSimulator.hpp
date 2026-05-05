// src/simulators/ParallelMultiAssetSimulator.hpp
#pragma once
#include <vector>
#include <omp.h>
#include "../stochasticProcess/MultiDimensionalProcess.hpp"
#include "../integrators/MultiDimensionalEulerMaruyama.hpp"
#include "../core/Random.hpp"

class ParallelMultiAssetSimulator {
public:
    // Returns results[path][asset] = terminal price
    static std::vector<std::vector<double>> simulate(
        const MultiDimensionalProcess& process,
        const std::vector<double>& x0,
        double T, double dt, int paths,
        int numThreads = 0)
    {
        int steps = static_cast<int>(T / dt);
        std::vector<std::vector<double>> results(paths);

        if (numThreads > 0)
            omp_set_num_threads(numThreads);

        int maxThreads = omp_get_max_threads();
        std::vector<Random> rngs;
        rngs.reserve(maxThreads);
        for (int i = 0; i < maxThreads; i++)
            rngs.emplace_back(42 + i * 1000);

        #pragma omp parallel for schedule(dynamic, 64)
        for (int i = 0; i < paths; i++) {
            int tid = omp_get_thread_num();
            Random& rng = rngs[tid];

            std::vector<double> x = x0;
            double t = 0.0;
            for (int j = 0; j < steps; j++) {
                x = MultiDimensionalEulerMaruyama::step(process, x, t, dt, rng);
                t += dt;
            }
            results[i] = x;
        }
        return results;
    }
};