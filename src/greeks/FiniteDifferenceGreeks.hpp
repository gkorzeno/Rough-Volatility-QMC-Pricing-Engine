// src/greeks/FiniteDifferenceGreeks.hpp
#pragma once
#include <vector>
#include <cmath>
#include <functional>
#include "../payoffs/Payoff.hpp"
#include "../core/Random.hpp"
#include "../core/OpenMPCompat.hpp"
#include "Greeks.hpp"

class FiniteDifferenceGreeks {
private:
    // Simulate GBM paths with given parameters, return terminal prices
    // Uses fixed seeds per path for common random numbers
    static std::vector<double> simulateGBM(
        double S0, double r, double sigma,
        double T, double dt, int paths, int baseSeed = 42)
    {
        int steps = static_cast<int>(T / dt);
        std::vector<double> results(paths);

        int maxThreads = omp_get_max_threads();
        std::vector<Random> rngs;
        for (int i = 0; i < maxThreads; i++)
            rngs.emplace_back(baseSeed + i * 1000);

        #pragma omp parallel for schedule(dynamic, 64)
        for (int i = 0; i < paths; i++) {
            Random& rng = rngs[omp_get_thread_num()];
            double S = S0, t = 0.0;
            for (int j = 0; j < steps; j++) {
                double z = rng.normal();
                S = S * std::exp((r - 0.5*sigma*sigma)*dt
                                 + sigma*std::sqrt(dt)*z);
                t += dt;
            }
            results[i] = S;
        }
        return results;
    }

    static double mcPrice(
        const std::vector<double>& paths,
        const Payoff& payoff,
        double r, double T)
    {
        double sum = 0.0;
        for (double S : paths) sum += payoff(S);
        return std::exp(-r * T) * sum / paths.size();
    }

public:
    static Greeks compute(
        const Payoff& payoff,
        double S0, double r, double sigma, double T,
        double dt, int paths,
        double dS    = 0.0,   // 0 = auto
        double dSig  = 0.0,
        double dR    = 0.0,
        double dT    = 0.0)
    {
        // Auto bump sizes if not specified
        if (dS   == 0.0) dS   = S0 * 0.01;       // 1% of S
        if (dSig == 0.0) dSig = sigma * 0.01;     // 1% of sigma
        if (dR   == 0.0) dR   = 0.0001;           // 1 basis point
        if (dT   == 0.0) dT   = 1.0 / 365.0;     // 1 day

        Greeks g;

        // Base price
        auto base = simulateGBM(S0, r, sigma, T, dt, paths);
        g.price   = mcPrice(base, payoff, r, T);

        // Delta and Gamma — bump S up and down
        auto up   = simulateGBM(S0 + dS, r, sigma, T, dt, paths);
        auto down = simulateGBM(S0 - dS, r, sigma, T, dt, paths);
        double pUp   = mcPrice(up,   payoff, r, T);
        double pDown = mcPrice(down, payoff, r, T);

        g.delta = (pUp - pDown) / (2.0 * dS);
        g.gamma = (pUp - 2.0 * g.price + pDown) / (dS * dS);

        // Vega — bump sigma
        auto vegaUp   = simulateGBM(S0, r, sigma + dSig, T, dt, paths);
        auto vegaDown = simulateGBM(S0, r, sigma - dSig, T, dt, paths);
        g.vega = (mcPrice(vegaUp, payoff, r, T)
                - mcPrice(vegaDown, payoff, r, T)) / (2.0 * dSig);

        // Rho — bump r
        // auto rhoUp   = simulateGBM(S0, r + dR, sigma, T, dt, paths);
        // auto rhoDown = simulateGBM(S0, r - dR, sigma, T, dt, paths);
        // g.rho = (mcPrice(rhoUp, payoff, r + dR, T)
        //        - mcPrice(rhoDown, payoff, r - dR, T)) / (2.0 * dR);

        // Replace rho block in FiniteDifferenceGreeks::compute
        double dR_bump = 0.01;  // 100 basis points — large enough to clear MC noise
        auto rhoUp   = simulateGBM(S0, r + dR_bump, sigma, T, dt, paths);
        auto rhoDown = simulateGBM(S0, r - dR_bump, sigma, T, dt, paths);
        g.rho = (mcPrice(rhoUp,   payoff, r + dR_bump, T)
            - mcPrice(rhoDown, payoff, r - dR_bump, T))
            / (2.0 * dR_bump);

        // Theta — reduce T by one day
        // auto thetaPath = simulateGBM(S0, r, sigma, T - dT, dt, paths);
        // g.theta = (mcPrice(thetaPath, payoff, r, T - dT) - g.price) / dT
        //           / 365.0;  // per calendar day

        // Theta is reported per day in Greeks.hpp and BlackScholes.
        auto thetaPath = simulateGBM(S0, r, sigma, T - dT, dt, paths);
        g.theta = (mcPrice(thetaPath, payoff, r, T - dT) - g.price) / dT / 365.0;

        return g;
    }
};
