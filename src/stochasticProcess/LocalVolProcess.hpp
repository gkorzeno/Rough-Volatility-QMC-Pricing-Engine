// src/stochasticProcess/LocalVolProcess.hpp
#pragma once
#include "StochasticProcess.hpp"
#include "../surface/LocalVolSurface.hpp"

// dS = r*S*dt + sigma_loc(S,t)*S*dW
class LocalVolProcess : public StochasticProcess {
private:
    const LocalVolSurface& lvSurface;
    double r;

public:
    LocalVolProcess(const LocalVolSurface& lv, double r_)
        : lvSurface(lv), r(r_) {}

    double drift(double x, double t) const override {
        return r * std::max(x, 1e-8);
    }

    double diffusion(double x, double t) const override {
        const double spot = std::max(x, 1e-8);
        return lvSurface.localVol(spot, t) * spot;
    }

    double diffusionDerivative(double x, double t) const override {
        return 0.0;  // numerical differentiation would be needed
    }
};
