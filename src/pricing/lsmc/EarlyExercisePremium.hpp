#pragma once
#include <vector>
#include <cmath>
#include <string>
#include "AmericanPricer.hpp"
#include "BermudanPricer.hpp"
#include "../../pricing/BlackScholes.hpp"

// Computes and decomposes the early exercise premium.
// Compares American/Bermudan prices against:
//   1. European MC price (same paths, no early exercise)
//   2. Black-Scholes European analytical price
//   3. Black-Scholes American approximation (Barone-Adesi-Whaley)
class EarlyExercisePremium {
public:
    struct Decomposition {
        double americanPrice;
        double europeanMC;
        double europeanBS;
        double americanBAW;       // Barone-Adesi-Whaley approximation
        double premiumVsEuroMC;   // American - European MC
        double premiumVsEuroBS;   // American - European BS
        double premiumVsBAW;      // American MC - American BAW (model impact)
        double roughVolImpact;    // American(rough) - American(flatVol)
    };

    // Compute the premium for a put under rough Bergomi vs flat vol
    static Decomposition analyzePut(
        const RoughBergomiModel::Parameters& roughParams,
        double spot, double K,
        double rate, double flatVol, double maturity,
        int nSteps, int nPaths,
        const AmericanPricer::Config& cfg = AmericanPricer::Config(),
        unsigned int seed = 42)
    {
        PutPayoff put(K);

        // American price under rough Bergomi
        auto roughResult = AmericanPricer::price(
            roughParams, put, spot, rate, maturity,
            nSteps, nPaths, cfg, seed);

        // American price under flat GBM (rBergomi with H→0.5, eta→0 approximation)
        // We approximate flat vol by setting xi0 = flatVol^2, eta ≈ 0
        RoughBergomiModel::Parameters flatParams = roughParams;
        flatParams.xi0   = flatVol * flatVol;
        flatParams.eta   = 1e-4;   // near-zero vol-of-vol → near-GBM
        flatParams.hurst = 0.499;  // near-Brownian

        auto flatResult = AmericanPricer::price(
            flatParams, put, spot, rate, maturity,
            nSteps, nPaths, cfg, seed + 1000);

        // European BS price
        const double euroBS = BlackScholes::putPrice(spot, K, rate, flatVol, maturity);

        // Barone-Adesi-Whaley approximation
        const double baw = baroneAdesiWhaleyPut(spot, K, rate, flatVol, maturity);

        return {
            roughResult.price,
            roughResult.europeanPrice,
            euroBS,
            baw,
            roughResult.earlyExercisePremium,
            roughResult.price - euroBS,
            roughResult.price - baw,
            roughResult.price - flatResult.price   // rough vol impact on EEP
        };
    }

private:
    struct PutPayoff : public Payoff {
        double K;
        PutPayoff(double k) : K(k) {}

        double operator()(double S) const override {
            return std::max(K - S, 0.0);
        }
    };

    // Barone-Adesi-Whaley approximation for American put
    static double baroneAdesiWhaleyPut(
        double S, double K, double r, double sigma, double T)
    {
        if (T <= 1e-8) return std::max(K - S, 0.0);
        const double eu = BlackScholes::putPrice(S, K, r, sigma, T);

        // Critical commodity price approximation
        const double M   = 2.0 * r / (sigma * sigma);
        const double N   = 2.0 * r / (sigma * sigma); // r_b=0
        const double q2  = (-(N - 1) - std::sqrt((N-1)*(N-1) + 4*M/
                            (1.0 - std::exp(-r*T)))) / 2.0;

        // Find critical spot S* via Newton (simplified)
        double Sstar = K * 0.5;
        for (int iter = 0; iter < 50; ++iter) {
            const double put_S = BlackScholes::putPrice(Sstar, K, r, sigma, T);
            const double d1    = (std::log(Sstar/K) + (r + 0.5*sigma*sigma)*T)
                               / (sigma * std::sqrt(T));
            const double N_d1  = 0.5 * std::erfc(d1 / std::sqrt(2.0));
            const double lhs   = K - Sstar;
            const double rhs   = put_S + (1.0 + 1.0/q2) * (-N_d1) * Sstar;
            if (std::abs(lhs - rhs) < 1e-6 * K) break;
            const double denom = -1.0 - (1.0/q2) * (-N_d1) + (1.0 + 1.0/q2) *
                                  std::exp(-0.5*d1*d1) / (sigma*std::sqrt(2.0*3.14159265358979323846*T));
            Sstar -= (lhs - rhs) / denom;
            Sstar = std::max(Sstar, 1e-4);
        }

        if (S <= Sstar)
            return K - S;  // exercise immediately

        const double put_star = BlackScholes::putPrice(Sstar, K, r, sigma, T);
        const double d1_star  = (std::log(Sstar/K) + (r+0.5*sigma*sigma)*T)
                              / (sigma*std::sqrt(T));
        const double N_d1_star = 0.5 * std::erfc(d1_star / std::sqrt(2.0));
        const double A2 = -(Sstar / q2) * (1.0 + N_d1_star);  // simplified

        return eu + A2 * std::pow(S / Sstar, q2);
    }
};