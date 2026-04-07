#pragma once
#include <cmath>
#include <vector>
#include "../core/Random.hpp"
#include "../core/OpenMPCompat.hpp"
#include "Greeks.hpp"

class AdjointGreeks {
public:
    // Reverse-mode sensitivities for a vanilla European call under GBM.
    // This first implementation targets first-order Greeks efficiently:
    // price, delta, vega, and rho.
    static Greeks computeCall(
        double K,
        double S0, double r, double sigma, double T,
        double dt, int paths)
    {
        const int steps = static_cast<int>(T / dt);
        Greeks g = {};

        double sumPrice = 0.0;
        double sumDelta = 0.0;
        double sumVega = 0.0;
        double sumRho = 0.0;

        const int maxThreads = omp_get_max_threads();
        std::vector<Random> rngs;
        rngs.reserve(maxThreads);
        for (int i = 0; i < maxThreads; ++i) {
            rngs.emplace_back(42 + i * 1000);
        }

        #pragma omp parallel for schedule(dynamic, 64) \
            reduction(+:sumPrice, sumDelta, sumVega, sumRho)
        for (int i = 0; i < paths; ++i) {
            Random& rng = rngs[omp_get_thread_num()];

            std::vector<double> spots(steps + 1, 0.0);
            std::vector<double> normals(steps, 0.0);
            spots[0] = S0;

            for (int j = 0; j < steps; ++j) {
                const double z = rng.normal();
                normals[j] = z;

                const double growth = std::exp((r - 0.5 * sigma * sigma) * dt
                                             + sigma * std::sqrt(dt) * z);
                spots[j + 1] = spots[j] * growth;
            }

            const double ST = spots[steps];
            const double intrinsic = std::max(ST - K, 0.0);
            const double disc = std::exp(-r * T);
            const double price = disc * intrinsic;
            sumPrice += price;

            if (intrinsic <= 0.0) {
                continue;
            }

            double barR = -T * price;
            double barSigma = 0.0;
            double barSpotNext = disc;

            for (int j = steps - 1; j >= 0; --j) {
                const double current = spots[j];
                const double next = spots[j + 1];

                const double dNext_dCurrent = next / current;
                const double dNext_dR = next * dt;
                const double dNext_dSigma =
                    next * (-sigma * dt + std::sqrt(dt) * normals[j]);

                barR += barSpotNext * dNext_dR;
                barSigma += barSpotNext * dNext_dSigma;
                barSpotNext *= dNext_dCurrent;
            }

            sumDelta += barSpotNext;
            sumVega += barSigma;
            sumRho += barR;
        }

        g.price = sumPrice / paths;
        g.delta = sumDelta / paths;
        g.vega = sumVega / paths;
        g.rho = sumRho / paths;
        g.gamma = 0.0;
        g.theta = 0.0;
        return g;
    }
};
