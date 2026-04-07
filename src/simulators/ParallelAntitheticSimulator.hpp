// src/simulators/ParallelAntitheticSimulator.hpp
#pragma once
#include <vector>
#include <stdexcept>
#include <omp.h>
#include "../stochasticProcess/StochasticProcess.hpp"
#include "../core/Random.hpp"

template<typename Integrator>
class ParallelAntitheticSimulator {
public:
    static std::vector<double> simulate(
        const StochasticProcess& process,
        double x0, double T, double dt, int paths,
        int numThreads = 0)
    {
        if (paths % 2 != 0)
            throw std::runtime_error("ParallelAntitheticSimulator requires even paths");

        int steps = static_cast<int>(T / dt);
        std::vector<double> results(paths);

        if (numThreads > 0)
            omp_set_num_threads(numThreads);

        int maxThreads = omp_get_max_threads();
        std::vector<Random> rngs;
        rngs.reserve(maxThreads);
        for (int i = 0; i < maxThreads; i++)
            rngs.emplace_back(42 + i * 1000);

        #pragma omp parallel for schedule(dynamic, 32)
        for (int i = 0; i < paths; i += 2) {
            int tid    = omp_get_thread_num();
            Random& rng = rngs[tid];

            double x1 = x0, x2 = x0, t = 0.0;
            for (int j = 0; j < steps; j++) {
                double z = rng.normal();
                x1 = Integrator::stepWithZ(process, x1, t, dt,  z);
                x2 = Integrator::stepWithZ(process, x2, t, dt, -z);
                t += dt;
            }
            results[i]   = x1;
            results[i+1] = x2;
        }

        return results;
    }
};