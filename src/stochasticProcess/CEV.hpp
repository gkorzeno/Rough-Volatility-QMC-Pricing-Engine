// src/stochasticProcess/CEV.hpp
#pragma once
#include "StochasticProcess.hpp"
#include <cmath>
#include <algorithm>

class CEV : public StochasticProcess {
private:
    double mu, sigma, beta;

public:
    CEV(double mu_, double sigma_, double beta_)
        : mu(mu_), sigma(sigma_), beta(beta_) {}

    double drift(double x, double t) const override {
        return mu * x;
    }

    double diffusion(double x, double t) const override {
        double xPos = std::max(x, 0.0);
        return sigma * std::pow(xPos, beta);
    }

    double diffusionDerivative(double x, double t) const override {
        double xPos = std::max(x, 1e-10);
        return sigma * beta * std::pow(xPos, beta - 1.0);
    }

    // CEV must use reflection — override floor
    double applyBoundary(double x) const {
        return std::max(x, 0.0);  // absorb at zero
    }

    bool useLogStep() const { return false; }
};