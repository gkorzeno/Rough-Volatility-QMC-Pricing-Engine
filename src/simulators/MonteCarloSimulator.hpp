#pragma once
#include <vector>
#include "../integrators/EulerMaruyama.hpp"

template<typename Integrator>
class MonteCarloSimulator {
public:
    static std::vector<double> simulate(
        const StochasticProcess& process,
        double x0,
        double T,
        double dt,
        int paths)
    {
        int steps = static_cast<int>(T / dt);
        std::vector<double> results(paths);

        Random rng;

        for (int i = 0; i < paths; i++) {
            double x = x0;
            double t = 0.0;

            for (int j = 0; j < steps; j++) {
                x = Integrator::step(process, x, t, dt, rng);
                // x = EulerMaruyama::step(process, x, t, dt, rng);
                //Double Jump:
                // if (process.hasJumps()) {
                //     int n = rng.poisson(process.jumpIntensity(x, t) * dt);

                //     for (int k = 0; k < n; k++) {
                //         double jump = process.sampleJumpSize(rng);
                //         x *= (1.0 + jump);
                //     }
                // }
                t += dt;
            }

            results[i] = x;
        }

        return results;
    }
};