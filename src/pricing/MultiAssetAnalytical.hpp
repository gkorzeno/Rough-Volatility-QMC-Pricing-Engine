// src/pricing/MultiAssetAnalytical.hpp
#pragma once
#include <cmath>
#include <vector>

class MultiAssetAnalytical {
    static double N(double x) {
        return 0.5 * std::erfc(-x / std::sqrt(2.0));
    }

public:
    // Margrabe formula: price of Max(S1 - S2, 0) at T
    // Assumes zero strike spread option with no dividends
    static double margrabePrice(
        double S1, double S2,
        double sigma1, double sigma2,
        double rho, double T)
    {
        double sigmaBar = std::sqrt(
            sigma1*sigma1 + sigma2*sigma2 - 2*rho*sigma1*sigma2);

        double d1 = (std::log(S1/S2) + 0.5*sigmaBar*sigmaBar*T)
                    / (sigmaBar * std::sqrt(T));
        double d2 = d1 - sigmaBar * std::sqrt(T);

        return S1 * N(d1) - S2 * N(d2);
    }

    // Approximate basket call via moment matching
    // Matches first two moments of weighted sum to lognormal
    static double basketCallApprox(
        const std::vector<double>& S0,
        const std::vector<double>& weights,
        const std::vector<double>& sigma,
        const std::vector<std::vector<double>>& rho,
        double r, double K, double T)
    {
        int n = S0.size();

        // First moment: E[basket] = sum(w_i * S_i * exp(r*T))
        double M1 = 0.0;
        for (int i = 0; i < n; i++)
            M1 += weights[i] * S0[i] * std::exp(r * T);

        // Second moment
        double M2 = 0.0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                M2 += weights[i] * weights[j]
                    * S0[i] * S0[j]
                    * std::exp((2*r + rho[i][j]*sigma[i]*sigma[j]) * T);

        // Match to lognormal
        double sigmaBasket = std::sqrt(std::log(M2 / (M1*M1)) / T);
        double S0basket    = M1 * std::exp(-r * T);

        // Black-Scholes with matched parameters
        double d1 = (std::log(S0basket/K) + (r + 0.5*sigmaBasket*sigmaBasket)*T)
                    / (sigmaBasket * std::sqrt(T));
        double d2 = d1 - sigmaBasket * std::sqrt(T);

        double Nd1 = N(d1), Nd2 = N(d2);
        return std::exp(-r*T) * (S0basket * std::exp(r*T) * Nd1 - K * Nd2);
    }
};