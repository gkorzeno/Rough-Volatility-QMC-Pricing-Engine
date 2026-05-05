#pragma once
#include "StochasticProcess.hpp"

class GeometricBrownianMotion : public StochasticProcess {
private:
    double mu;
    double sigma;

public:
    GeometricBrownianMotion(double mu_, double sigma_)
        : mu(mu_), sigma(sigma_) {}

    double drift(double x, double t) const override {
        return mu * x;
    }

    double diffusion(double x, double t) const override {
        return sigma * x;
    }

    // GeometricBrownianMotion.hpp
    double diffusionDerivative(double x, double t) const override {
        return sigma;  // d/dx(sigma*x) = sigma
    }

};