#pragma once
#include <functional>
#include <stdexcept>
#include <vector>

class FractionalSDE {
public:
    struct Path {
        std::vector<double> times;
        std::vector<double> states;
        std::vector<double> drivers;
    };

    static Path simulateEuler(
        double x0,
        double T,
        const std::vector<double>& fractionalDriverPath,
        const std::function<double(double, double)>& drift,
        const std::function<double(double, double)>& diffusion)
    {
        const int steps = static_cast<int>(fractionalDriverPath.size());
        if (steps == 0 || T <= 0.0)
            throw std::runtime_error("FractionalSDE: invalid horizon/driver");

        const double dt = T / steps;
        Path path;
        path.times.resize(steps + 1, 0.0);
        path.states.resize(steps + 1, 0.0);
        path.drivers.resize(steps + 1, 0.0);
        path.states[0] = x0;

        double prevDriver = 0.0;
        for (int i = 0; i < steps; ++i) {
            const double t = i * dt;
            const double dBH = fractionalDriverPath[i] - prevDriver;
            prevDriver = fractionalDriverPath[i];

            path.times[i + 1] = (i + 1) * dt;
            path.drivers[i + 1] = fractionalDriverPath[i];
            path.states[i + 1] =
                path.states[i]
                + drift(path.states[i], t) * dt
                + diffusion(path.states[i], t) * dBH;
        }

        return path;
    }
};
