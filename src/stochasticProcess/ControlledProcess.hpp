// src/stochasticProcess/ControlledProcess.hpp
#pragma once
#include "StochasticProcess.hpp"
#include "../control/Policy.hpp"

class ControlledProcess : public StochasticProcess {
public:
    // Controlled versions — subclasses implement these
    virtual double drift(double x, double u, double t)     const = 0;
    virtual double diffusion(double x, double u, double t) const = 0;

    // Satisfy the base interface by evaluating with u=0
    // (uncontrolled fallback so existing simulator still works)
    double drift(double x, double t)     const override { return drift(x, 0.0, t); }
    double diffusion(double x, double t) const override { return diffusion(x, 0.0, t); }
};