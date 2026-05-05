// src/stochasticProcess/SABR.hpp
#pragma once
#include "CholeskyProcess.hpp"
#include <cmath>
#include <algorithm>

// State vector: x = [F, alpha]
// F     = forward price
// alpha = stochastic volatility
class SABR : public CholeskyProcess {
private:
    double beta;   // elasticity (0=normal, 1=lognormal)
    double nu;     // vol of vol
    double rho;    // correlation

public:
    SABR(double beta_, double nu_, double rho_)
        : beta(beta_), nu(nu_), rho(rho_) {}

    std::vector<double> drift(
        const std::vector<double>& x, double t) const override
    {
        return {0.0, 0.0};  // SABR is driftless under forward measure
    }

    // Returns L directly (Cholesky factored) since we manage correlation
    std::vector<std::vector<double>> diffusion(
        const std::vector<double>& x, double t) const override
    {
        double F     = std::max(x[0], 1e-8);
        double alpha = std::max(x[1], 1e-8);
        double Fb    = std::pow(F, beta);

        return {
            {alpha * Fb,                              0.0                           },
            {rho * nu * alpha, nu * alpha * std::sqrt(1.0 - rho*rho)}
        };
    }

    // Hagan et al. implied volatility approximation — for validation
    double impliedVol(double F, double K, double alpha,
                      double T, double r) const
    {
        if (std::abs(F - K) < 1e-10) {  // ATM formula
            double FK_mid = std::pow(F*K, (1-beta)/2.0);
            double term1  = 1.0
                + ((1-beta)*(1-beta)/24.0 * alpha*alpha / std::pow(F, 2-2*beta)
                +  0.25 * rho*beta*nu*alpha / std::pow(F, 1-beta)
                +  (2-3*rho*rho)/24.0 * nu*nu) * T;
            return alpha / FK_mid * term1;
        }

        double logFK   = std::log(F/K);
        double FK_mid  = std::pow(F*K, (1-beta)/2.0);
        double z       = nu/alpha * FK_mid * logFK;
        double x_z     = std::log((std::sqrt(1 - 2*rho*z + z*z) + z - rho)
                                  / (1 - rho));
        double term1   = alpha / (FK_mid * (1 + (1-beta)*(1-beta)/24.0
                         * logFK*logFK
                         + std::pow(1-beta, 4)/1920.0 * logFK*logFK*logFK*logFK));
        double term2   = (std::abs(x_z) < 1e-10) ? 1.0 : z / x_z;
        double term3   = 1.0
            + ((1-beta)*(1-beta)/24.0 * alpha*alpha / std::pow(F*K, 1-beta)
            +  0.25 * rho*beta*nu*alpha / std::pow(F*K, (1-beta)/2.0)
            +  (2-3*rho*rho)/24.0 * nu*nu) * T;

        return term1 * term2 * term3;
    }
};