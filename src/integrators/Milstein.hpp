// src/integrators/Milstein.hpp
#pragma once
#include <cmath>
#include "../stochasticProcess/StochasticProcess.hpp"
#include "../core/Random.hpp"

class Milstein {
public:
    static double step(
        const StochasticProcess& process,
        double x,
        double t,
        double dt,
        Random& rng)
    {
        double mu      = process.drift(x, t);
        double sig     = process.diffusion(x, t);
        double sig_    = process.diffusionDerivative(x, t);
        double z       = rng.normal();
        double sqrtDt  = std::sqrt(dt);

        return x
            + mu * dt
            + sig * sqrtDt * z
            + 0.5 * sig * sig_ * (z*z - 1.0) * dt;  // Milstein correction
    }
};