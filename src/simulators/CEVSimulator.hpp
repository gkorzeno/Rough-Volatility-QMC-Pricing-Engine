// src/simulators/CEVSimulator.hpp
#pragma once
#include <vector>
#include <omp.h>
#include "../stochasticProcess/CEV.hpp"
#include "../integrators/CEVEulerMaruyama.hpp"
#include "../core/Random.hpp"

class CEVSimulator {
public:
    static std::vector<double> simulate(
        const CEV& process,
        double x0, double T, double dt, int paths)
    {
        int steps = static_cast<int>(T / dt);
        std::vector<double> results(paths);

        int maxThreads = omp_get_max_threads();
        std::vector<Random> rngs;
        for (int i = 0; i < maxThreads; i++)
            rngs.emplace_back(42 + i * 1000);

        #pragma omp parallel for schedule(dynamic, 64)
        for (int i = 0; i < paths; i++) {
            Random& rng = rngs[omp_get_thread_num()];
            double x = x0, t = 0.0;
            for (int j = 0; j < steps; j++) {
                x = CEVEulerMaruyama::step(process, x, t, dt, rng);
                t += dt;
            }
            results[i] = x;
        }
        return results;
    }
};