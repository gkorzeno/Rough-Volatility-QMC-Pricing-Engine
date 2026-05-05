// src/calibration/ImpliedVolSolver.hpp
#pragma once
#include <cmath>
#include <stdexcept>
#include "../pricing/BlackScholes.hpp"

class ImpliedVolSolver {
public:
    // Newton-Raphson inversion of BS call price
    static double solve(
        double marketPrice,
        double S, double K, double r, double T,
        double sigmaInit = 0.2,
        int    maxIter   = 100,
        double tol       = 1e-8)
    {
        // Bounds check — price must be within arbitrage bounds
        double intrinsic = std::max(S - K * std::exp(-r * T), 0.0);
        if (marketPrice <= intrinsic)
            throw std::runtime_error("Price below intrinsic — no implied vol exists");
        if (marketPrice >= S)
            throw std::runtime_error("Price above spot — no implied vol exists");

        double sigma = sigmaInit;

        for (int i = 0; i < maxIter; i++) {
            double price = BlackScholes::callPrice(S, K, r, sigma, T);
            double vega  = BlackScholes::vega(S, K, r, sigma, T);

            double diff  = price - marketPrice;
            if (std::abs(diff) < tol) return sigma;
            if (std::abs(vega) < 1e-12) break;  // flat vega — switch to bisection

            sigma -= diff / vega;
            sigma  = std::max(sigma, 1e-6);  // keep positive
        }

        // Fallback: bisection if Newton diverges
        return bisection(marketPrice, S, K, r, T);
    }

private:
    static double bisection(
        double target,
        double S, double K, double r, double T,
        double lo = 1e-6, double hi = 5.0,
        double tol = 1e-8)
    {
        for (int i = 0; i < 200; i++) {
            double mid   = 0.5 * (lo + hi);
            double price = BlackScholes::callPrice(S, K, r, mid, T);
            if (std::abs(price - target) < tol) return mid;
            (price < target ? lo : hi) = mid;
        }
        return 0.5 * (lo + hi);
    }
};