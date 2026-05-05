#pragma once
#include "../core/Random.hpp"

class StochasticProcess {
public:
    //Drift Term
    virtual double drift(double x, double t) const = 0;
    //Diffusion Term
    virtual double diffusion(double x, double t) const = 0;

    //Jump Diffusion Extension
    virtual bool hasJumps() const { return false; }
    virtual double jumpIntensity(double x, double t) const { return 0.0; }
    virtual double sampleJumpSize(Random& rng) const { return 0.0; }
    virtual double diffusionDerivative(double x, double t) const { return 0.0; }

    virtual ~StochasticProcess() {}
};