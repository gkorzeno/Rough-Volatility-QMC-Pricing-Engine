// src/simulators/ControlledMonteCarloSimulator.hpp
#pragma once
#include <vector>
#include "../integrators/ControlledEulerMaruyama.hpp"

class ControlledMonteCarloSimulator {
public:
    static std::vector<double> simulate(
        const ControlledProcess& process,
        const Policy& policy,
        double x0,
        double T,
        double dt,
        int paths)
    {
        int steps = static_cast<int>(T / dt);
        std::vector<double> results(paths);
        Random rng;

        for (int i = 0; i < paths; i++) {
            double x = x0, t = 0.0;
            for (int j = 0; j < steps; j++) {
                x = ControlledEulerMaruyama::step(process, policy, x, t, dt, rng);
                t += dt;
            }
            results[i] = x;
        }
        return results;
    }
};