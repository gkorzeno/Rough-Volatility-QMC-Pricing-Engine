// src/simulators/PathSimulator.hpp — stores entire paths, not just terminal values
#pragma once
#include <vector>
#include "../integrators/EulerMaruyama.hpp"

template<typename Integrator>
class PathSimulator {
public:
    // Returns paths[i][j] = price of path i at step j
    static std::vector<std::vector<double>> simulate(
        const StochasticProcess& process,
        double x0, double T, double dt, int paths)
    {
        int steps = static_cast<int>(T / dt);
        std::vector<std::vector<double>> allPaths(paths,
            std::vector<double>(steps + 1));
        Random rng;

        for (int i = 0; i < paths; i++) {
            allPaths[i][0] = x0;
            double x = x0, t = 0.0;
            for (int j = 0; j < steps; j++) {
                x = Integrator::step(process, x, t, dt, rng);
                t += dt;
                allPaths[i][j+1] = x;
            }
        }
        return allPaths;
    }
};