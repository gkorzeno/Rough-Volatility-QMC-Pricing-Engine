// src/integrators/ControlledEulerMaruyama.hpp
#pragma once
#include <cmath>
#include "../stochasticProcess/ControlledProcess.hpp"
#include "../control/Policy.hpp"
#include "../core/Random.hpp"

class ControlledEulerMaruyama {
public:
    static double step(
        const ControlledProcess& process,
        const Policy& policy,
        double x,
        double t,
        double dt,
        Random& rng)
    {
        double u         = policy.control(x, t);
        double drift     = process.drift(x, u, t);
        double diffusion = process.diffusion(x, u, t);
        double z         = rng.normal();

        return x + drift * dt + diffusion * std::sqrt(dt) * z;
    }
};