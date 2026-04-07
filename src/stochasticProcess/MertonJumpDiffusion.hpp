// MertonJumpDiffusion.hpp
#pragma once
#include "StochasticProcess.hpp"
#include "../core/Random.hpp"

class MertonJumpDiffusion : public StochasticProcess {
private:
    double mu;       // drift
    double sigma;    // diffusion volatility
    double lambda;   // jump arrival rate (jumps per year)
    double muJ;      // mean log jump size
    double sigmaJ;   // std dev of log jump size
    double kappa;

public:
    MertonJumpDiffusion(double mu_, double sigma_,
                        double lambda_, double muJ_, double sigmaJ_)
        : mu(mu_), sigma(sigma_),
          lambda(lambda_), muJ(muJ_), sigmaJ(sigmaJ_), 
          kappa(std::exp(muJ_ + 0.5 * sigmaJ_ * sigmaJ_) - 1.0) {}

    // Drift is compensated: subtract expected jump contribution
    // so the process remains a martingale (risk-neutral pricing)
    double drift(double x, double t) const override {
        return (mu - lambda * kappa) * x;
    }

    double diffusion(double x, double t) const override {
        return sigma * x;
    }

    bool hasJumps() const override { return true; }

    double jumpIntensity(double x, double t) const override {
        return lambda;
    }

    // Returns the multiplicative jump factor (e^J - 1), so caller does x *= (1 + factor)
    double sampleJumpSize(Random& rng) const override {
        double J = muJ + sigmaJ * rng.normal();
        return std::exp(J) - 1.0;
    }

    double theoreticalMean(double x0, double T) const {
        // E[X_T] = x0 * exp(mu * T)
        // The compensated drift ensures this holds exactly
        return x0 * std::exp(mu * T);
    }

    double theoreticalVariance(double x0, double T) const {
        // Var[X_T] = x0^2 * exp(2*mu*T) * (exp(v*T) - 1)
        // where v = sigma^2 + lambda*(exp(2*muJ + sigmaJ^2)*(exp(sigmaJ^2) - 1))
        double kappa2 = std::exp(2*muJ + sigmaJ*sigmaJ) * (std::exp(sigmaJ*sigmaJ) - 1.0);
        double v = sigma*sigma + lambda * kappa2;
        return x0*x0 * std::exp(2*mu*T) * (std::exp(v*T) - 1.0);
    }

    // // MertonJumpDiffusion.hpp
    // double diffusionDerivative(double x, double t) const override {
    //     return sigma;  // same as GBM
    // }
};