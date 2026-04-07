// src/greeks/PathwiseGreeks.hpp
#pragma once
#include <vector>
#include <cmath>
#include "../core/Random.hpp"
#include "../core/OpenMPCompat.hpp"
#include "Greeks.hpp"

class PathwiseGreeks {
public:
    // Returns delta and vega via pathwise estimators
    // More efficient than finite difference for smooth payoffs
    static Greeks compute(
        double K,
        double S0, double r, double sigma, double T,
        double dt, int paths)
    {
        int steps = static_cast<int>(T / dt);
        Greeks g  = {};

        double sumPrice = 0.0, sumDelta = 0.0;
        double sumVega  = 0.0, sumRho   = 0.0;

        int maxThreads = omp_get_max_threads();
        std::vector<Random> rngs;
        for (int i = 0; i < maxThreads; i++)
            rngs.emplace_back(42 + i * 1000);

        // Use reduction for thread safety
        #pragma omp parallel for schedule(dynamic,64) \
            reduction(+:sumPrice, sumDelta, sumVega, sumRho)
        for (int i = 0; i < paths; i++) {
            Random& rng = rngs[omp_get_thread_num()];

            double S         = S0;
            double sumZ      = 0.0;  // accumulated noise for vega
            double t         = 0.0;

            for (int j = 0; j < steps; j++) {
                double z  = rng.normal();
                sumZ     += z;
                S = S * std::exp((r - 0.5*sigma*sigma)*dt
                                 + sigma*std::sqrt(dt)*z);
                t += dt;
            }

            if (S > K) {
                double disc = std::exp(-r * T);

                // Price
                sumPrice += disc * (S - K);

                // Delta: dPayoff/dS0 = S_T/S0 * 1_{S_T > K}
                sumDelta += disc * S / S0;

                // Vega: dPayoff/dsigma = S_T * (sumZ*sqrt(dt) - sigma*T)
                sumVega  += disc * S * (sumZ * std::sqrt(dt) - sigma * T);

                // Rho: dPayoff/dr = S_T * T
                sumRho   += disc * S * T;
            }
        }

        double disc = std::exp(-r * T);
        g.price = sumPrice / paths;
        g.delta = sumDelta / paths;
        g.vega  = sumVega  / paths;
        g.rho   = disc * sumRho / paths;

        // Gamma not available via pathwise for vanilla call
        // (payoff not twice differentiable at K)
        g.gamma = 0.0;
        g.theta = 0.0;

        return g;
    }
};
