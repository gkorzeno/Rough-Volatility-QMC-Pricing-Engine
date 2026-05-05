// src/stochasticProcess/Heston.hpp
#pragma once
#include "CholeskyProcess.hpp" //MultidimensionalProcess
#include <cmath>
#include <algorithm>

// State vector: x = [S, v]
class Heston : public CholeskyProcess { //MultidimensionalProcess
private:
    double mu;     // asset drift
    double kappa;  // mean reversion speed of variance
    double vbar;   // long-run variance
    double xi;     // volatility of variance ("vol of vol")
    double rho;    // correlation between S and v

public:
    Heston(double mu_, double kappa_, double vbar_,
           double xi_, double rho_)
        : mu(mu_), kappa(kappa_), vbar(vbar_),
          xi(xi_), rho(rho_) {}

    // x[0] = S (asset price)
    // x[1] = v (variance)
    std::vector<double> drift(
        const std::vector<double>& x, double t) const override
    {
        double S = x[0];
        double v = std::max(x[1], 0.0);  // floor variance

        return {
            mu * S,
            kappa * (vbar - v)
        };
    }

    // Covariance matrix of the noise:
    // [v*S^2,        rho*xi*v*S]
    // [rho*xi*v*S,   xi^2*v    ]
    std::vector<std::vector<double>> diffusion(
        const std::vector<double>& x, double t) const override
    {
        double S   = x[0];
        double v   = std::max(x[1], 0.0);
        double sqv = std::sqrt(v);

        return {
            {sqv * S,        0.0        },
            {rho * xi * sqv, xi * sqv * std::sqrt(1 - rho*rho)}
        };
    }

    // Feller condition: if 2*kappa*vbar > xi^2, variance stays positive
    bool fellerConditionSatisfied() const {
        return 2 * kappa * vbar > xi * xi;
    }

    // Theoretical mean of log(S_T) under risk-neutral measure
    double theoreticalLogMean(double S0, double v0, double T) const {
        return std::log(S0) + (mu - 0.5 * vbar) * T;  // approximate
    }
};