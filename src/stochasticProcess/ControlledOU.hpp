// src/stochasticProcess/ControlledOU.hpp
#pragma once
#include "ControlledProcess.hpp"

// dX = (theta*(mu - x) + u) dt + sigma dW
// Control u adds a direct push — optimal u minimizes E[integral of (x^2 + r*u^2) dt]
class ControlledOU : public ControlledProcess {
private:
    double theta, mu, sigma;

public:
    ControlledOU(double theta_, double mu_, double sigma_)
        : theta(theta_), mu(mu_), sigma(sigma_) {}

    double drift(double x, double u, double t) const override {
        return theta * (mu - x) + u;
    }

    double diffusion(double x, double u, double t) const override {
        return sigma;
    }
};