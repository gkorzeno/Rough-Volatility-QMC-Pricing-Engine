// src/integrators/CEVEulerMaruyama.hpp
#pragma once
#include <cmath>
#include "../stochasticProcess/CEV.hpp"
#include "../core/Random.hpp"

class CEVEulerMaruyama {
public:
    static double step(
        const CEV& process,
        double x,
        double t,
        double dt,
        Random& rng)
    {
        double z    = rng.normal();
        double mu   = process.drift(x, t);
        double sig  = process.diffusion(x, t);

        double next_x = x + mu * dt + sig * std::sqrt(dt) * z;

        // Reflect at zero — Andersen absorption scheme
        // If next_x < 0, use full truncation
        return std::max(next_x, 0.0);
    }

    static double stepWithZ(
        const CEV& process,
        double x,
        double t,
        double dt,
        double z)
    {
        double mu  = process.drift(x, t);
        double sig = process.diffusion(x, t);
        return std::max(x + mu * dt + sig * std::sqrt(dt) * z, 0.0);
    }
};