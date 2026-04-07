#pragma once
#include <vector>
#include "../integrators/MultiDimensionalEulerMaruyama.hpp"

class MultiDimensionalMonteCarloSimulator {
public:

    static std::vector<std::vector<double>> simulate(
        const MultiDimensionalProcess& process,
        const std::vector<double>& x0,
        double T,
        double dt,
        int paths)
    {
        int steps = static_cast<int>(T / dt);
        std::vector<std::vector<double>> results(paths);

        Random rng;

        for (int i = 0; i < paths; i++) {

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