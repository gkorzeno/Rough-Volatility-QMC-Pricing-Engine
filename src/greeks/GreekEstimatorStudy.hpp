#pragma once
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
#include "../core/BrownianBridge.hpp"
#include "../core/Random.hpp"
#include "../core/SobolSequence.hpp"
#include "../pricing/BlackScholes.hpp"

class GreekEstimatorStudy {
public:
    enum Sampler {
        MC,
        QMC,
        RQMC
    };

    struct RunGreeks {
        Greeks pathwise = {};
        Greeks likelihoodRatio = {};
        Greeks aad = {};
        double pwDeltaStdErr = 0.0;
        double pwVegaStdErr = 0.0;
        double lrDeltaStdErr = 0.0;
        double lrVegaStdErr = 0.0;
        double aadDeltaStdErr = 0.0;
        double aadVegaStdErr = 0.0;
    };

    struct EstimatorSummary {
        double deltaMean = 0.0;
        double deltaRmse = 0.0;
        double deltaStability = std::numeric_limits<double>::quiet_NaN();
        double deltaInnerStdErr = 0.0;
        double vegaMean = 0.0;
        double vegaRmse = 0.0;
        double vegaStability = std::numeric_limits<double>::quiet_NaN();
        double vegaInnerStdErr = 0.0;
    };

    struct Summary {
        EstimatorSummary pathwise;
        EstimatorSummary likelihoodRatio;
        EstimatorSummary aad;
    };

    static RunGreeks runSingle(
        Sampler sampler,
        double K,
        double S0,
        double r,
        double sigma,
        double T,
        double dt,
        int paths,
        std::uint64_t seed = 42,
        const std::string& directionFile = "docs/new-joe-kuo-6.21201.txt",
        bool useBrownianBridge = true,
        bool useOwenScrambling = false,
        int rqmcReplicates = 8)
    {
        const int steps = static_cast<int>(T / dt + 1e-12);
        const double sqrtDt = std::sqrt(dt);
        const double disc = std::exp(-r * T);

        Accumulator pw;
        Accumulator lr;
        Accumulator aad;

        auto consumePath = [&](const std::vector<double>& z) {
            std::vector<double> spots(steps + 1, 0.0);
            spots[0] = S0;

            double sumZ = 0.0;
            for (int j = 0; j < steps; ++j) {
                sumZ += z[j];
                const double growth = std::exp((r - 0.5 * sigma * sigma) * dt
                                             + sigma * sqrtDt * z[j]);
                spots[j + 1] = spots[j] * growth;
            }

            const double ST = spots[steps];
            const double intrinsic = std::max(ST - K, 0.0);
            const double priceContribution = disc * intrinsic;
            const double WT = sumZ * sqrtDt;
            const double ZT = WT / std::sqrt(T);

            double pwDelta = 0.0;
            double pwVega = 0.0;
            if (intrinsic > 0.0) {
                pwDelta = disc * ST / S0;
                pwVega = disc * ST * (WT - sigma * T);
            }

            const double scoreDelta = WT / (S0 * sigma * T);
            const double scoreVega = ((ZT * ZT - 1.0) / sigma) - std::sqrt(T) * ZT;
            const double lrDelta = priceContribution * scoreDelta;
            const double lrVega = priceContribution * scoreVega;

            double aadDelta = 0.0;
            double aadVega = 0.0;
            if (intrinsic > 0.0) {
                double barSigma = 0.0;
                double barSpotNext = disc;
                for (int j = steps - 1; j >= 0; --j) {
                    const double current = spots[j];
                    const double next = spots[j + 1];
                    const double dNext_dCurrent = next / current;
                    const double dNext_dSigma = next * (-sigma * dt + sqrtDt * z[j]);
                    barSigma += barSpotNext * dNext_dSigma;
                    barSpotNext *= dNext_dCurrent;
                }
                aadDelta = barSpotNext;
                aadVega = barSigma;
            }

            pw.add(priceContribution, pwDelta, pwVega);
            lr.add(priceContribution, lrDelta, lrVega);
            aad.add(priceContribution, aadDelta, aadVega);
        };

        if (sampler == MC) {
            Random rng(static_cast<unsigned int>(seed & 0xFFFFFFFFULL));
            std::vector<double> z(steps, 0.0);
            for (int i = 0; i < paths; ++i) {
                for (int j = 0; j < steps; ++j)
                    z[j] = rng.normal();
                consumePath(z);
            }
        } else if (sampler == QMC) {
            SobolSequence sobol(steps, directionFile);
            std::vector<double> z(steps, 0.0);
            for (int i = 0; i < paths; ++i) {
                auto u = sobol.next();
                fillNormalsFromUniforms(z, u, dt, useBrownianBridge);
                consumePath(z);
            }
        } else {
            if (rqmcReplicates <= 1 || paths % rqmcReplicates != 0)
                throw std::runtime_error("GreekEstimatorStudy: RQMC paths must be divisible by replicates and replicates > 1");

            const int pointsPerReplicate = paths / rqmcReplicates;
            std::vector<double> shift(steps, 0.0);
            std::vector<double> z(steps, 0.0);
            for (int rep = 0; rep < rqmcReplicates; ++rep) {
                const std::uint64_t repSeed = seed + static_cast<std::uint64_t>(rep) * 0x9E3779B97F4A7C15ULL;
                SobolSequence sobol = useOwenScrambling
                    ? SobolSequence(steps, directionFile, repSeed, 12)
                    : SobolSequence(steps, directionFile);

                if (!useOwenScrambling) {
                    std::uint64_t shiftSeed = splitmix64(repSeed);
                    for (int j = 0; j < steps; ++j) {
                        shiftSeed = splitmix64(shiftSeed);
                        constexpr double inv2_53 = 1.0 / 9007199254740992.0;
                        shift[j] = static_cast<double>((shiftSeed >> 11) & ((1ULL << 53) - 1)) * inv2_53;
                    }
                }

                for (int i = 0; i < pointsPerReplicate; ++i) {
                    auto u = sobol.next();
                    if (!useOwenScrambling) {
                        for (int j = 0; j < steps; ++j) {
                            u[j] += shift[j];
                            if (u[j] >= 1.0) u[j] -= 1.0;
                        }
                    }
                    fillNormalsFromUniforms(z, u, dt, useBrownianBridge);
                    consumePath(z);
                }
            }
        }

        RunGreeks result;
        result.pathwise = pw.mean(paths);
        result.likelihoodRatio = lr.mean(paths);
        result.aad = aad.mean(paths);
        result.pwDeltaStdErr = pw.stdErrDelta(paths);
        result.pwVegaStdErr = pw.stdErrVega(paths);
        result.lrDeltaStdErr = lr.stdErrDelta(paths);
        result.lrVegaStdErr = lr.stdErrVega(paths);
        result.aadDeltaStdErr = aad.stdErrDelta(paths);
        result.aadVegaStdErr = aad.stdErrVega(paths);
        return result;
    }

    static Summary summarize(
        Sampler sampler,
        double K,
        double S0,
        double r,
        double sigma,
        double T,
        double dt,
        int paths,
        int outerRuns,
        const std::string& directionFile = "docs/new-joe-kuo-6.21201.txt",
        bool useBrownianBridge = true,
        bool useOwenScrambling = false,
        int rqmcReplicates = 8,
        std::uint64_t baseSeed = 42)
    {
        const Greeks bs = BlackScholes::greeks(S0, K, r, sigma, T);

        std::vector<double> pwDelta, pwVega, lrDelta, lrVega, aadDelta, aadVega;
        std::vector<double> pwDeltaSe, pwVegaSe, lrDeltaSe, lrVegaSe, aadDeltaSe, aadVegaSe;

        const int runs = (sampler == QMC) ? 1 : outerRuns;
        for (int run = 0; run < runs; ++run) {
            const std::uint64_t seed = baseSeed + static_cast<std::uint64_t>(run) * 7919ULL;
            const auto estimate = runSingle(
                sampler, K, S0, r, sigma, T, dt, paths, seed,
                directionFile, useBrownianBridge, useOwenScrambling, rqmcReplicates);

            pwDelta.push_back(estimate.pathwise.delta);
            pwVega.push_back(estimate.pathwise.vega);
            lrDelta.push_back(estimate.likelihoodRatio.delta);
            lrVega.push_back(estimate.likelihoodRatio.vega);
            aadDelta.push_back(estimate.aad.delta);
            aadVega.push_back(estimate.aad.vega);

            pwDeltaSe.push_back(estimate.pwDeltaStdErr);
            pwVegaSe.push_back(estimate.pwVegaStdErr);
            lrDeltaSe.push_back(estimate.lrDeltaStdErr);
            lrVegaSe.push_back(estimate.lrVegaStdErr);
            aadDeltaSe.push_back(estimate.aadDeltaStdErr);
            aadVegaSe.push_back(estimate.aadVegaStdErr);
        }

        Summary summary;
        summary.pathwise = summarizeEstimator(pwDelta, pwVega, pwDeltaSe, pwVegaSe, bs.delta, bs.vega);
        summary.likelihoodRatio = summarizeEstimator(lrDelta, lrVega, lrDeltaSe, lrVegaSe, bs.delta, bs.vega);
        summary.aad = summarizeEstimator(aadDelta, aadVega, aadDeltaSe, aadVegaSe, bs.delta, bs.vega);
        return summary;
    }

    static std::string samplerLabel(Sampler sampler, bool useOwenScrambling) {
        if (sampler == MC) return "MC";
        if (sampler == QMC) return "QMC";
        return useOwenScrambling ? "RQMC (Owen)" : "RQMC (Shift)";
    }

private:
    struct Accumulator {
        double sumPrice = 0.0;
        double sumDelta = 0.0;
        double sumDeltaSq = 0.0;
        double sumVega = 0.0;
        double sumVegaSq = 0.0;

        void add(double price, double delta, double vega) {
            sumPrice += price;
            sumDelta += delta;
            sumDeltaSq += delta * delta;
            sumVega += vega;
            sumVegaSq += vega * vega;
        }

        Greeks mean(int n) const {
            Greeks g = {};
            g.price = sumPrice / n;
            g.delta = sumDelta / n;
            g.vega = sumVega / n;
            return g;
        }

        double stdErrDelta(int n) const {
            return stdErr(sumDelta, sumDeltaSq, n);
        }

        double stdErrVega(int n) const {
            return stdErr(sumVega, sumVegaSq, n);
        }

        static double stdErr(double sum, double sumSq, int n) {
            const double mean = sum / n;
            const double var = std::max(0.0, (sumSq / n) - mean * mean);
            return std::sqrt(var / n);
        }
    };

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

        const double pLow = 0.02425;
        const double pHigh = 1.0 - pLow;
        double q = 0.0, r = 0.0;
        if (u < pLow) {
            q = std::sqrt(-2.0 * std::log(u));
            return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
                   ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
        }
        if (u <= pHigh) {
            q = u - 0.5;
            r = q * q;
            return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
                   (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
        }
        q = std::sqrt(-2.0 * std::log(1.0 - u));
        return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
                ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
    }

    static void fillNormalsFromUniforms(std::vector<double>& z, const std::vector<double>& u, double dt, bool useBrownianBridge) {
        for (size_t j = 0; j < z.size(); ++j)
            z[j] = normalInverse(u[j]);
        if (useBrownianBridge) {
            const auto dW = BrownianBridge::build(z, dt);
            const double invSqrtDt = 1.0 / std::sqrt(dt);
            for (size_t j = 0; j < z.size(); ++j)
                z[j] = dW[j] * invSqrtDt;
        }
    }

    static std::uint64_t splitmix64(std::uint64_t x) {
        x += 0x9E3779B97F4A7C15ULL;
        x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
        x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
        return x ^ (x >> 31);
    }

    static EstimatorSummary summarizeEstimator(
        const std::vector<double>& delta,
        const std::vector<double>& vega,
        const std::vector<double>& deltaSe,
        const std::vector<double>& vegaSe,
        double targetDelta,
        double targetVega)
    {
        EstimatorSummary s;
        s.deltaMean = mean(delta);
        s.deltaRmse = rmse(delta, targetDelta);
        s.deltaStability = stability(delta);
        s.deltaInnerStdErr = mean(deltaSe);
        s.vegaMean = mean(vega);
        s.vegaRmse = rmse(vega, targetVega);
        s.vegaStability = stability(vega);
        s.vegaInnerStdErr = mean(vegaSe);
        return s;
    }

    static double mean(const std::vector<double>& x) {
        return std::accumulate(x.begin(), x.end(), 0.0) / x.size();
    }

    static double rmse(const std::vector<double>& x, double target) {
        double sq = 0.0;
        for (double v : x) {
            const double d = v - target;
            sq += d * d;
        }
        return std::sqrt(sq / x.size());
    }

    static double stability(const std::vector<double>& x) {
        if (x.size() < 2)
            return std::numeric_limits<double>::quiet_NaN();
        const double m = mean(x);
        double sq = 0.0;
        for (double v : x) {
            const double d = v - m;
            sq += d * d;
        }
        return std::sqrt(sq / (x.size() - 1));
    }
};
