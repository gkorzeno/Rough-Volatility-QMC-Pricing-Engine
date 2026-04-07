// src/greeks/LikelihoodRatioGreeks.hpp
#pragma once
#include <vector>
#include <cmath>
#include "../core/Random.hpp"
#include "../core/OpenMPCompat.hpp"
#include "../payoffs/Payoff.hpp"
#include "Greeks.hpp"

class LikelihoodRatioGreeks {
public:
    static Greeks compute(
        const Payoff& payoff,
        double S0, double r, double sigma, double T,
        double dt, int paths)
    {
        int steps = static_cast<int>(T / dt);
        Greeks g  = {};

        double sumPrice = 0.0, sumDelta = 0.0;
        double sumGamma = 0.0, sumVega  = 0.0;

        int maxThreads = omp_get_max_threads();
        std::vector<Random> rngs;
        for (int i = 0; i < maxThreads; i++)
            rngs.emplace_back(42 + i * 1000);

        #pragma omp parallel for schedule(dynamic, 64) \
            reduction(+:sumPrice, sumDelta, sumGamma, sumVega)
        for (int i = 0; i < paths; i++) {
            Random& rng = rngs[omp_get_thread_num()];

            double S = S0;
            double t = 0.0;

            for (int j = 0; j < steps; j++) {
                double z = rng.normal();
                S = S * std::exp((r - 0.5*sigma*sigma)*dt
                                 + sigma*std::sqrt(dt)*z);
                t += dt;
            }

            double payout = payoff(S);
            double disc   = std::exp(-r * T);

            // Compute W_T from terminal S rather than accumulating z
            // W_T = (log(S_T/S0) - (r - 0.5*sigma^2)*T) / sigma
            double logRet = std::log(S / S0);
            double mu_T   = (r - 0.5*sigma*sigma) * T;
            double W_T    = (logRet - mu_T) / sigma;
            double Z      = W_T / std::sqrt(T);

            // Score for delta: W_T / (S0 * sigma * T)
            double scoreDelta = W_T / (S0 * sigma * T);

            // LR second-order weight:
            // d2 log p / dS0^2 + (d log p / dS0)^2
            double scoreGamma = (Z*Z / (sigma*sigma*T)
                               - Z / (sigma * std::sqrt(T))
                               - 1.0) / (S0 * S0);

            // Score for vega in lognormal model
            double scoreVega = (Z*Z - 1.0) / sigma - std::sqrt(T) * Z;

            sumPrice += disc * payout;
            sumDelta += disc * payout * scoreDelta;
            sumGamma += disc * payout * scoreGamma;
            sumVega  += disc * payout * scoreVega;
        }

        g.price = sumPrice / paths;
        g.delta = sumDelta / paths;
        g.gamma = sumGamma / paths;
        g.vega  = sumVega  / paths;
        g.theta = 0.0;
        g.rho   = 0.0;

        return g;
    }
};
