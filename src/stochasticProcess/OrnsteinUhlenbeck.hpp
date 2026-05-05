#pragma once
#include "StochasticProcess.hpp"
#include <cmath>

class OrnsteinUhlenbeck : public StochasticProcess {
private:
    double theta;
    double mu;
    double sigma;

public:
    OrnsteinUhlenbeck(double theta_, double mu_, double sigma_)
        : theta(theta_), mu(mu_), sigma(sigma_) {}

    double drift(double x, double t) const override {
        return theta * (mu - x);
    }

    double diffusion(double x, double t) const override {
        return sigma;
    }

    double theoreticalVariance(double t) const {
        return (sigma * sigma) / (2 * theta) * (1 - std::exp(-2 * theta * t));
    }

    double theoreticalMean(double x0, double t) const {
        return mu + (x0 - mu) * std::exp(-theta * t);
    }

    // OrnsteinUhlenbeck.hpp
    // double diffusionDerivative(double x, double t) const override {
    //     return 0.0;  // diffusion is constant
    // }
};