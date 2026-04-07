#pragma once
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>
#include <omp.h>
#include "../stochasticProcess/StochasticProcess.hpp"
#include "../payoffs/Payoff.hpp"
#include "../core/SobolSequence.hpp"
#include "../core/BrownianBridge.hpp"

template<typename Integrator>
class ParallelRQMCSimulator {
public:
    struct Result {
        double price;
        double stdErr;      // standard error across randomized replicates
        int replicates;
        int pointsPerReplicate;
    };

    static Result priceEuropean(
        const StochasticProcess& process,
        const Payoff& payoff,
        double x0, double r, double T, double dt,
        int pointsPerReplicate,
        int replicates,
        const std::string& directionFile,
        bool useBrownianBridge = true,
        bool useOwenScrambling = false,
        int owenDepth = 12,
        std::uint64_t baseSeed = 123456789ULL)
    {
        if (pointsPerReplicate <= 0 || replicates <= 1)
            throw std::runtime_error("ParallelRQMC: pointsPerReplicate>0 and replicates>1 required");
        if (dt <= 0.0 || T <= 0.0)
            throw std::runtime_error("ParallelRQMC: T and dt must be positive");

        int steps = static_cast<int>(T / dt);
        if (steps <= 0)
            throw std::runtime_error("ParallelRQMC: invalid time grid");

        std::vector<double> replicateMeans(replicates, 0.0);
        double sqrtDt = std::sqrt(dt);
        double disc   = std::exp(-r * T);

        #pragma omp parallel for schedule(dynamic, 1)
        for (int rep = 0; rep < replicates; ++rep) {
            std::uint64_t repSeed = baseSeed + static_cast<std::uint64_t>(rep) * 0x9E3779B97F4A7C15ULL;
            SobolSequence sobol = useOwenScrambling
                ? SobolSequence(steps, directionFile, repSeed, owenDepth)
                : SobolSequence(steps, directionFile);

            // If not using Owen scrambling, add a randomized digital shift
            // (standard RQMC). If using Owen scrambling, the scramble itself
            // already randomizes, so we skip the digital shift.
            std::vector<double> shift;
            if (!useOwenScrambling) {
                shift.resize(steps);
                std::uint64_t seed = splitmix64(baseSeed + static_cast<std::uint64_t>(rep) * 0x9E3779B97F4A7C15ULL);
                for (int j = 0; j < steps; ++j) {
                    seed = splitmix64(seed);
                    constexpr double inv2_53 = 1.0 / 9007199254740992.0; // 2^53
                    shift[j] = static_cast<double>((seed >> 11) & ((1ULL << 53) - 1)) * inv2_53;
                }
            }

            double sumPayoff = 0.0;
            std::vector<double> z(steps), incOrZ(steps);

            for (int n = 0; n < pointsPerReplicate; ++n) {
                auto u = sobol.next();

                for (int j = 0; j < steps; ++j) {
                    double uj = u[j];
                    if (!useOwenScrambling) {
                        uj += shift[j];
                        if (uj >= 1.0) uj -= 1.0;
                    }
                    z[j] = normalInverse(uj);
                }

                if (useBrownianBridge) {
                    auto dW = BrownianBridge::build(z, dt);
                    for (int j = 0; j < steps; ++j)
                        incOrZ[j] = dW[j] / sqrtDt;
                } else {
                    incOrZ = z;
                }

                double x = x0;
                double t = 0.0;
                for (int j = 0; j < steps; ++j) {
                    x = Integrator::stepWithZ(process, x, t, dt, incOrZ[j]);
                    t += dt;
                }

                sumPayoff += payoff(x);
            }

            replicateMeans[rep] = disc * (sumPayoff / pointsPerReplicate);
        }

        double mean = 0.0;
        for (double m : replicateMeans) mean += m;
        mean /= replicates;

        double sq = 0.0;
        for (double m : replicateMeans) {
            double d = m - mean;
            sq += d * d;
        }
        double sampleVar = sq / (replicates - 1);
        double se = std::sqrt(sampleVar / replicates);

        return {mean, se, replicates, pointsPerReplicate};
    }

private:
    static std::uint64_t splitmix64(std::uint64_t x) {
        x += 0x9E3779B97F4A7C15ULL;
        x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
        x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
        return x ^ (x >> 31);
    }

    // Acklam-style inverse CDF approximation.
    static double normalInverse(double u) {
        if (u <= 0.0) u = 1e-15;
        if (u >= 1.0) u = 1.0 - 1e-15;

        static const double a[] = {
            -3.969683028665376e+01,  2.209460984245205e+02,
            -2.759285104469687e+02,  1.383577518672690e+02,
            -3.066479806614716e+01,  2.506628277459239e+00
        };
        static const double b[] = {
            -5.447609879822406e+01,  1.615858368580409e+02,
            -1.556989798598866e+02,  6.680131188771972e+01,
            -1.328068155288572e+01
        };
        static const double c[] = {
            -7.784894002430293e-03, -3.223964580411365e-01,
            -2.400758277161838e+00, -2.549732539343734e+00,
             4.374664141464968e+00,  2.938163982698783e+00
        };
        static const double d[] = {
             7.784695709041462e-03,  3.224671290700398e-01,
             2.445134137142996e+00,  3.754408661907416e+00
        };

        const double pLow  = 0.02425;
        const double pHigh = 1.0 - pLow;

        double q, r;
        if (u < pLow) {
            q = std::sqrt(-2.0 * std::log(u));
            return (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
                   ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0);
        }
        if (u <= pHigh) {
            q = u - 0.5;
            r = q * q;
            return (((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5]) * q /
                   (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1.0);
        }

        q = std::sqrt(-2.0 * std::log(1.0 - u));
        return -(((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
                ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0);
    }
};
