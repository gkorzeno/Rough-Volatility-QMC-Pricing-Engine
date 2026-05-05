// src/simulators/MultiDimPathSimulator.hpp
#pragma once
#include <vector>
#include "../stochasticProcess/MultiDimensionalProcess.hpp"
#include "../integrators/MultiDimensionalEulerMaruyama.hpp"
#include "../core/Random.hpp"

// Returns paths[path][step][dim]
class MultiDimPathSimulator {
public:
    static std::vector<std::vector<std::vector<double>>> simulate(
        const MultiDimensionalProcess& process,
        const std::vector<double>& x0,
        double T, double dt, int paths)
    {
        int steps = static_cast<int>(T / dt);
        int dim   = x0.size();

        std::vector<std::vector<std::vector<double>>> allPaths(
            paths, std::vector<std::vector<double>>(
                steps + 1, std::vector<double>(dim)));

        Random rng;

        for (int i = 0; i < paths; i++) {
            allPaths[i][0] = x0;
            std::vector<double> x = x0;
            double t = 0.0;

            for (int j = 0; j < steps; j++) {
                x = MultiDimensionalEulerMaruyama::step(process, x, t, dt, rng);
                allPaths[i][j+1] = x;
                t += dt;
            }
        }
        return allPaths;
    }
};