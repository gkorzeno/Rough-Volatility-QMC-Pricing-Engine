// TwoDimOU.hpp
#pragma once
#include "MultiDimensionalProcess.hpp"
#include <cmath>

class TwoDimensionalOU : public MultiDimensionalProcess {
private:
    double theta, mu1, mu2, sigma1, sigma2, rho;

public:
    TwoDimensionalOU(double theta_, double mu1_, double mu2_,
             double sigma1_, double sigma2_, double rho_)
        : theta(theta_), mu1(mu1_), mu2(mu2_),
          sigma1(sigma1_), sigma2(sigma2_), rho(rho_) {}

    std::vector<double> drift(
        const std::vector<double>& x, double t) const override
    {
        return {
            -theta * (x[0] - mu1),
            -theta * (x[1] - mu2)
        };
    }

    // Returns covariance matrix — Cholesky handles the rest
    std::vector<std::vector<double>> diffusion(
        const std::vector<double>& x, double t) const override
    {
        return {
            {sigma1 * sigma1,        rho * sigma1 * sigma2},
            {rho * sigma1 * sigma2,  sigma2 * sigma2      }
        };
    }

    // Theoretical mean: mu + (x0 - mu) * exp(-theta * t)
    std::vector<double> theoreticalMean(
        const std::vector<double>& x0, double t) const override
    {
        return {
            mu1 + (x0[0] - mu1) * std::exp(-theta * t),
            mu2 + (x0[1] - mu2) * std::exp(-theta * t)
        };
    }

    // Theoretical marginal variance: sigma^2 / (2*theta) * (1 - exp(-2*theta*t))
    std::vector<double> theoreticalVariance(double t) const override
    {
        double v1 = (sigma1*sigma1) / (2*theta) * (1 - std::exp(-2*theta*t));
        double v2 = (sigma2*sigma2) / (2*theta) * (1 - std::exp(-2*theta*t));
        return {v1, v2};
    }
};