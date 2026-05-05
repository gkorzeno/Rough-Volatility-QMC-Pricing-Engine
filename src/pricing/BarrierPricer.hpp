// src/pricing/BarrierPricer.hpp
#pragma once
#include <vector>
#include <cmath>
#include "../stochasticProcess/GeometricBrownianMotion.hpp"
#include "../barriers/Barrier.hpp"
#include "../barriers/BarrierCrossing.hpp"
#include "../payoffs/Payoff.hpp"
#include "../core/Random.hpp"
#include <omp.h>

class BarrierPricer {
public:
    struct Result {
        double price;
        double stdErr;
        double fractionKnockedOut;
    };

    static Result price(
        const GeometricBrownianMotion& process,
        const Payoff& payoff,
        const Barrier& barrier,
        double S0, double r, double sigma,
        double T, double dt, int paths,
        bool useBridgeDetection = true)
    {
        int steps = static_cast<int>(T / dt);
        std::vector<double> payoffs(paths, 0.0);
        std::vector<int>    knockedOut(paths, 0);

        int maxThreads = omp_get_max_threads();
        std::vector<Random> rngs;
        rngs.reserve(maxThreads);
        for (int i = 0; i < maxThreads; i++)
            rngs.emplace_back(42 + i * 1000);

        #pragma omp parallel for schedule(dynamic, 64)
        for (int i = 0; i < paths; i++) {
            int tid     = omp_get_thread_num();
            Random& rng = rngs[tid];

            double S       = S0;
            double t       = 0.0;
            bool   alive   = true;
            bool   touched = false;

            for (int j = 0; j < steps && alive; j++) {
                double S_prev = S;
                double z      = rng.normal();
                S = S_prev * std::exp(
                    (r - 0.5*sigma*sigma)*dt + sigma*std::sqrt(dt)*z);

                // ── Barrier detection ──────────────────────────────────────
                bool crossed = false;
                if (useBridgeDetection) {
                    bool s0above = (S_prev >= barrier.H);
                    bool s1above = (S      >= barrier.H);

                    if (s0above != s1above) {
                        // Endpoints on opposite sides — definite cross
                        crossed = true;
                    } else if (barrier.isUp() && !s0above) {
                        // Up barrier: both below — bridge check
                        double prob = BarrierCrossing::crossingProbability(
                            S_prev, S, barrier.H, sigma, dt);
                        // crossed = rng.uniform() < prob;
                        // if (prob > 0.0 && rng.uniform() < prob)
                        //     crossed = true;
                        
                        if (prob < 1e-12) {
                            crossed = false;   // skip tiny probabilities (numerical noise)
                        } else {
                            crossed = rng.uniform() < prob;
                        }
                    } else if (!barrier.isUp() && s0above) {
                        // Down barrier: both above — bridge check
                        double prob = BarrierCrossing::crossingProbability(
                            S_prev, S, barrier.H, sigma, dt);
                        // crossed = rng.uniform() < prob;
                        // if (prob > 0.0 && rng.uniform() < prob)
                        //     crossed = true;

                        if (prob < 1e-12) {
                            crossed = false;   // skip tiny probabilities (numerical noise)
                        } else {
                            crossed = rng.uniform() < prob;
                        }
                    }
                    // Both above up-barrier or both below down-barrier:
                    // no crossing possible, crossed stays false
                } else {
                    crossed = barrier.isUp()
                        ? (S_prev < barrier.H && S >= barrier.H)
                        : (S_prev > barrier.H && S <= barrier.H);
                }
                // ──────────────────────────────────────────────────────────

                if (crossed) {
                    touched = true;
                    if (barrier.isOut()) {
                        alive = false;
                        knockedOut[i] = 1;
                        break;  // ← stop simulating immediately after knock-out
                    }
                }

                t += dt;
            }

            if (barrier.isOut()) {
                payoffs[i] = alive ? payoff(S) : barrier.rebate;
            } else {
                payoffs[i] = touched ? payoff(S) : 0.0;
            }
        }

        int    n    = paths;
        double mean = 0.0;
        for (double p : payoffs) mean += p;
        mean /= n;

        double sq_sum = 0.0;
        for (double p : payoffs) sq_sum += (p - mean) * (p - mean);
        double stdErr = std::sqrt(sq_sum / n) / std::sqrt(n);

        double knockFrac = 0.0;
        for (int k : knockedOut) knockFrac += k;
        knockFrac /= n;

        return {std::exp(-r*T) * mean, std::exp(-r*T) * stdErr, knockFrac};
    }
};