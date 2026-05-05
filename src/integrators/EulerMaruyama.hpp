#pragma once
#include <cmath>
#include "../stochasticProcess/StochasticProcess.hpp"
#include "../core/Random.hpp"

class EulerMaruyama {
public:
    static double step(
        const StochasticProcess& process,
        double x,
        double t,
        double dt,
        Random& rng)
    {
        double drift = process.drift(x, t);
        double diffusion = process.diffusion(x, t);
        double z = rng.normal();

        double next_x = x
            + drift * dt
            + diffusion * std::sqrt(dt) * z;

        //Jump Diffusion
        //Poisson number of jumps
        if (process.hasJumps()) {
            double intensity = process.jumpIntensity(x, t);
            int numJumps = rng.poisson(intensity * dt);  // see note below
            for (int k = 0; k < numJumps; k++) {
                next_x += next_x * process.sampleJumpSize(rng);
            }
        }

        return next_x;
    }

    static double stepWithZ(
        const StochasticProcess& process,
        double x,
        double t,
        double dt,
        double z)
    {
        double mu  = process.drift(x, t);
        double sig = process.diffusion(x, t);

        return x + mu * dt + sig * std::sqrt(dt) * z;
    }
};