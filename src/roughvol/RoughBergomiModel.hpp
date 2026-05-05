#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>
#include "../core/SobolSequence.hpp"
#include "../core/Random.hpp"
#include "../calibration/ImpliedVolSolver.hpp"
#include "../surface/VolSurface.hpp"

class RoughBergomiModel {
public:
    enum SamplingMethod {
        MC,
        QMC,
        RQMC
    };

    struct Parameters {
        double xi0;
        double eta;
        double rho;
        double hurst;
    };

    struct SurfaceSpec {
        double spot;
        double rate;
        double maturity;
        int steps;
        int paths;
        int rqmcReplicates;
        std::vector<double> strikes;
        std::vector<double> expiries;
    };

    struct PricingResult {
        double price;
        double stdErr;
        int paths;
        SamplingMethod method;
    };

    struct PathResult {
        std::vector<double> times;
        std::vector<double> spot;
        std::vector<double> variance;
        std::vector<double> volterra;
        std::vector<double> volBrownianIncrements;
        std::vector<double> orthogonalBrownianIncrements;
        std::vector<double> spotBrownianIncrements;
    };

    struct DriverCorrelationResult {
        double targetRho;
        double empiricalRho;
        int samples;
    };

    static std::vector<double> simulateVolterraProcess(
        double maturity,
        int steps,
        double hurst,
        Random& rng)
    {
        if (hurst <= 0.0 || hurst >= 0.5)
            throw std::runtime_error("RoughBergomiModel: H should be in (0, 0.5) for rough volatility");

        const double dt = maturity / steps;
        std::vector<double> dw(steps, 0.0);
        std::vector<double> y(steps, 0.0);

        for (int i = 0; i < steps; ++i)
            dw[i] = std::sqrt(dt) * rng.normal();

        for (int i = 0; i < steps; ++i) {
            double acc = 0.0;
            const double ti = (i + 1) * dt;
            for (int j = 0; j <= i; ++j) {
                const double left = j * dt;
                const double right = (j + 1) * dt;
                const double mid = 0.5 * (left + right);
                const double lag = std::max(ti - mid, 1e-12);
                acc += std::pow(lag, hurst - 0.5) * dw[j];
            }
            y[i] = std::sqrt(2.0 * hurst) * acc;
        }

        return y;
    }

    static PathResult simulatePath(
        const Parameters& params,
        double spot,
        double rate,
        double maturity,
        int steps,
        Random& rng)
    {
        if (params.xi0 <= 0.0 || params.eta <= 0.0)
            throw std::runtime_error("RoughBergomiModel: xi0 and eta must be positive");
        if (std::abs(params.rho) >= 1.0)
            throw std::runtime_error("RoughBergomiModel: rho must lie in (-1,1)");

        const double dt = maturity / steps;
        const std::vector<double> y = simulateVolterraProcess(maturity, steps, params.hurst, rng);

        PathResult result;
        result.times.resize(steps + 1, 0.0);
        result.spot.resize(steps + 1, 0.0);
        result.variance.resize(steps + 1, params.xi0);
        result.volterra.resize(steps + 1, 0.0);
        result.volBrownianIncrements.resize(steps, 0.0);
        result.orthogonalBrownianIncrements.resize(steps, 0.0);
        result.spotBrownianIncrements.resize(steps, 0.0);
        result.spot[0] = spot;

        for (int i = 0; i < steps; ++i) {
            const double tNext = (i + 1) * dt;
            const double variance =
                params.xi0
                * std::exp(params.eta * y[i]
                         - 0.5 * params.eta * params.eta * std::pow(tNext, 2.0 * params.hurst));

            const double dW1 = std::sqrt(dt) * rng.normal();
            const double dW2 = std::sqrt(dt) * rng.normal();
            const double dB = params.rho * dW1 + std::sqrt(1.0 - params.rho * params.rho) * dW2;

            result.times[i + 1] = tNext;
            result.volterra[i + 1] = y[i];
            result.variance[i + 1] = variance;
            result.volBrownianIncrements[i] = dW1;
            result.orthogonalBrownianIncrements[i] = dW2;
            result.spotBrownianIncrements[i] = dB;
            result.spot[i + 1] =
                result.spot[i]
                * std::exp((rate - 0.5 * variance) * dt + std::sqrt(std::max(variance, 0.0)) * dB);
        }

        return result;
    }

    static PathResult simulatePathFromNormals(
        const Parameters& params,
        double spot,
        double rate,
        double maturity,
        int steps,
        const std::vector<double>& normals)
    {
        if (static_cast<int>(normals.size()) != 3 * steps)
            throw std::runtime_error("RoughBergomiModel: expected 3*steps normals");

        if (params.xi0 <= 0.0 || params.eta <= 0.0)
            throw std::runtime_error("RoughBergomiModel: xi0 and eta must be positive");
        if (std::abs(params.rho) >= 1.0)
            throw std::runtime_error("RoughBergomiModel: rho must lie in (-1,1)");

        const double dt = maturity / steps;
        std::vector<double> dw(steps, 0.0);
        for (int i = 0; i < steps; ++i)
            dw[i] = std::sqrt(dt) * normals[i];

        std::vector<double> y(steps, 0.0);
        for (int i = 0; i < steps; ++i) {
            double acc = 0.0;
            const double ti = (i + 1) * dt;
            for (int j = 0; j <= i; ++j) {
                const double left = j * dt;
                const double right = (j + 1) * dt;
                const double mid = 0.5 * (left + right);
                const double lag = std::max(ti - mid, 1e-12);
                acc += std::pow(lag, params.hurst - 0.5) * dw[j];
            }
            y[i] = std::sqrt(2.0 * params.hurst) * acc;
        }

        PathResult result;
        result.times.resize(steps + 1, 0.0);
        result.spot.resize(steps + 1, 0.0);
        result.variance.resize(steps + 1, params.xi0);
        result.volterra.resize(steps + 1, 0.0);
        result.volBrownianIncrements.resize(steps, 0.0);
        result.orthogonalBrownianIncrements.resize(steps, 0.0);
        result.spotBrownianIncrements.resize(steps, 0.0);
        result.spot[0] = spot;

        for (int i = 0; i < steps; ++i) {
            const double tNext = (i + 1) * dt;
            const double variance =
                params.xi0
                * std::exp(params.eta * y[i]
                         - 0.5 * params.eta * params.eta * std::pow(tNext, 2.0 * params.hurst));

            const double dW1 = std::sqrt(dt) * normals[steps + 2 * i];
            const double dW2 = std::sqrt(dt) * normals[steps + 2 * i + 1];
            const double dB = params.rho * dW1 + std::sqrt(1.0 - params.rho * params.rho) * dW2;

            result.times[i + 1] = tNext;
            result.volterra[i + 1] = y[i];
            result.variance[i + 1] = variance;
            result.volBrownianIncrements[i] = dW1;
            result.orthogonalBrownianIncrements[i] = dW2;
            result.spotBrownianIncrements[i] = dB;
            result.spot[i + 1] =
                result.spot[i]
                * std::exp((rate - 0.5 * variance) * dt + std::sqrt(std::max(variance, 0.0)) * dB);
        }
        return result;
    }

    static PricingResult priceEuropeanCall(
        const Parameters& params,
        double spot,
        double rate,
        double maturity,
        int steps,
        int paths,
        double strike,
        SamplingMethod method = MC,
        const std::string& directionFile = "docs/new-joe-kuo-6.21201.txt",
        bool useOwenScrambling = true,
        int rqmcReplicates = 8,
        std::uint64_t baseSeed = 42ULL)
    {
        if (paths <= 0)
            throw std::runtime_error("RoughBergomiModel: paths must be positive");

        const double disc = std::exp(-rate * maturity);
        const int dim = 3 * steps;

        if (method == MC) {
            Random rng(static_cast<unsigned int>(baseSeed));
            std::vector<double> payoffs;
            payoffs.reserve(paths);
            std::vector<double> normals(dim, 0.0);

            for (int i = 0; i < paths; ++i) {
                for (int j = 0; j < dim; ++j)
                    normals[j] = rng.normal();
                PathResult path = simulatePathFromNormals(params, spot, rate, maturity, steps, normals);
                payoffs.push_back(disc * std::max(path.spot.back() - strike, 0.0));
            }
            return summarizePayoffs(payoffs, method);
        }

        if (method == QMC) {
            SobolSequence sobol(dim, directionFile);
            std::vector<double> payoffs;
            payoffs.reserve(paths);
            std::vector<double> normals(dim, 0.0);
            for (int i = 0; i < paths; ++i) {
                const std::vector<double> u = sobol.next();
                for (int j = 0; j < dim; ++j)
                    normals[j] = normalInverse(u[j]);
                PathResult path = simulatePathFromNormals(params, spot, rate, maturity, steps, normals);
                payoffs.push_back(disc * std::max(path.spot.back() - strike, 0.0));
            }
            PricingResult result = summarizePayoffs(payoffs, method);
            result.stdErr = std::numeric_limits<double>::quiet_NaN();
            return result;
        }

        if (rqmcReplicates <= 1 || paths % rqmcReplicates != 0)
            throw std::runtime_error("RoughBergomiModel: RQMC paths must be divisible by replicates and replicates>1");

        const int pointsPerRep = paths / rqmcReplicates;
        std::vector<double> replicateMeans(rqmcReplicates, 0.0);
        std::vector<double> shift(dim, 0.0);
        std::vector<double> normals(dim, 0.0);

        for (int rep = 0; rep < rqmcReplicates; ++rep) {
            const std::uint64_t repSeed = baseSeed + static_cast<std::uint64_t>(rep) * 0x9E3779B97F4A7C15ULL;
            SobolSequence sobol = useOwenScrambling
                ? SobolSequence(dim, directionFile, repSeed, 12)
                : SobolSequence(dim, directionFile);

            if (!useOwenScrambling) {
                std::uint64_t shiftSeed = splitmix64(repSeed);
                for (int j = 0; j < dim; ++j) {
                    shiftSeed = splitmix64(shiftSeed);
                    constexpr double inv2_53 = 1.0 / 9007199254740992.0;
                    shift[j] = static_cast<double>((shiftSeed >> 11) & ((1ULL << 53) - 1)) * inv2_53;
                }
            }

            double payoffSum = 0.0;
            for (int n = 0; n < pointsPerRep; ++n) {
                std::vector<double> u = sobol.next();
                for (int j = 0; j < dim; ++j) {
                    double uj = u[j];
                    if (!useOwenScrambling) {
                        uj += shift[j];
                        if (uj >= 1.0) uj -= 1.0;
                    }
                    normals[j] = normalInverse(uj);
                }
                PathResult path = simulatePathFromNormals(params, spot, rate, maturity, steps, normals);
                payoffSum += disc * std::max(path.spot.back() - strike, 0.0);
            }
            replicateMeans[rep] = payoffSum / pointsPerRep;
        }

        PricingResult result;
        result.price = std::accumulate(replicateMeans.begin(), replicateMeans.end(), 0.0) / rqmcReplicates;
        double sq = 0.0;
        for (int i = 0; i < rqmcReplicates; ++i) {
            const double d = replicateMeans[i] - result.price;
            sq += d * d;
        }
        result.stdErr = std::sqrt(sq / std::max(1, rqmcReplicates - 1)) / std::sqrt(static_cast<double>(rqmcReplicates));
        result.paths = paths;
        result.method = method;
        return result;
    }

    static VolSurface generateVolSurface(
        const Parameters& params,
        const SurfaceSpec& spec,
        int seed = 42,
        SamplingMethod method = MC,
        bool useOwenScrambling = true)
    {
        if (spec.expiries.empty() || spec.strikes.empty())
            throw std::runtime_error("RoughBergomiModel: surface grid cannot be empty");
        const std::string directionFile = "docs/new-joe-kuo-6.21201.txt";
        std::vector<std::vector<double>> implied(
            spec.expiries.size(),
            std::vector<double>(spec.strikes.size(), 0.0));

        if (method == MC) {
            std::vector<std::vector<double>> priceSums(
                spec.expiries.size(),
                std::vector<double>(spec.strikes.size(), 0.0));

            Random rng(seed);

            for (int path = 0; path < spec.paths; ++path) {
                PathResult simulated = simulatePath(
                    params, spec.spot, spec.rate, spec.maturity, spec.steps, rng);

                for (std::size_t ei = 0; ei < spec.expiries.size(); ++ei) {
                    const double expiry = spec.expiries[ei];
                    const int idx = std::min(
                        static_cast<int>(std::ceil(expiry / spec.maturity * spec.steps)),
                        spec.steps);
                    const double terminal = simulated.spot[idx];
                    const double disc = std::exp(-spec.rate * expiry);

                    for (std::size_t ki = 0; ki < spec.strikes.size(); ++ki)
                        priceSums[ei][ki] += disc * std::max(terminal - spec.strikes[ki], 0.0);
                }
            }

            for (std::size_t ei = 0; ei < spec.expiries.size(); ++ei) {
                for (std::size_t ki = 0; ki < spec.strikes.size(); ++ki) {
                    implied[ei][ki] = clampAndSolve(
                        priceSums[ei][ki] / spec.paths,
                        spec,
                        ei,
                        ki,
                        params);
                }
            }
            return VolSurface(spec.strikes, spec.expiries, implied);
        }

        for (std::size_t ei = 0; ei < spec.expiries.size(); ++ei) {
            const double expiry = spec.expiries[ei];
            const double disc = std::exp(-spec.rate * expiry);
            const int localSteps = std::max(
                1,
                static_cast<int>(std::ceil(expiry / spec.maturity * spec.steps)));
            const int dim = 3 * localSteps;
            std::vector<double> priceSums(spec.strikes.size(), 0.0);

            if (method == QMC) {
                SobolSequence sobol(dim, directionFile);
                std::vector<double> normals(dim, 0.0);

                for (int n = 0; n < spec.paths; ++n) {
                    const std::vector<double> u = sobol.next();
                    for (int j = 0; j < dim; ++j)
                        normals[j] = normalInverse(u[j]);

                    const PathResult path = simulatePathFromNormals(
                        params, spec.spot, spec.rate, expiry, localSteps, normals);
                    const double terminal = path.spot.back();
                    for (std::size_t ki = 0; ki < spec.strikes.size(); ++ki)
                        priceSums[ki] += disc * std::max(terminal - spec.strikes[ki], 0.0);
                }
            } else {
                if (spec.rqmcReplicates <= 1 || spec.paths % spec.rqmcReplicates != 0)
                    throw std::runtime_error("RQMC: paths must be divisible by rqmcReplicates");

                const int pointsPerRep = spec.paths / spec.rqmcReplicates;
                std::vector<std::vector<double>> replicateMeans(
                    spec.rqmcReplicates,
                    std::vector<double>(spec.strikes.size(), 0.0));

                for (int rep = 0; rep < spec.rqmcReplicates; ++rep) {
                    const std::uint64_t repSeed =
                        static_cast<std::uint64_t>(seed)
                        + static_cast<std::uint64_t>(rep) * 0x9E3779B97F4A7C15ULL;

                    SobolSequence sobol = useOwenScrambling
                        ? SobolSequence(dim, directionFile, repSeed, 12)
                        : SobolSequence(dim, directionFile);

                    std::vector<double> shift(dim, 0.0);
                    if (!useOwenScrambling) {
                        std::uint64_t shiftSeed = splitmix64(repSeed);
                        constexpr double inv2_53 = 1.0 / 9007199254740992.0;
                        for (int j = 0; j < dim; ++j) {
                            shiftSeed = splitmix64(shiftSeed);
                            shift[j] = static_cast<double>(
                                (shiftSeed >> 11) & ((1ULL << 53) - 1)) * inv2_53;
                        }
                    }

                    std::vector<double> normals(dim, 0.0);
                    std::vector<double> strikeSum(spec.strikes.size(), 0.0);
                    for (int n = 0; n < pointsPerRep; ++n) {
                        std::vector<double> u = sobol.next();
                        for (int j = 0; j < dim; ++j) {
                            double uj = u[j];
                            if (!useOwenScrambling) {
                                uj += shift[j];
                                if (uj >= 1.0)
                                    uj -= 1.0;
                            }
                            normals[j] = normalInverse(uj);
                        }

                        const PathResult path = simulatePathFromNormals(
                            params, spec.spot, spec.rate, expiry, localSteps, normals);
                        const double terminal = path.spot.back();
                        for (std::size_t ki = 0; ki < spec.strikes.size(); ++ki)
                            strikeSum[ki] += disc * std::max(terminal - spec.strikes[ki], 0.0);
                    }

                    for (std::size_t ki = 0; ki < spec.strikes.size(); ++ki)
                        replicateMeans[rep][ki] = strikeSum[ki] / pointsPerRep;
                }

                for (int rep = 0; rep < spec.rqmcReplicates; ++rep)
                    for (std::size_t ki = 0; ki < spec.strikes.size(); ++ki)
                        priceSums[ki] += replicateMeans[rep][ki];

                for (std::size_t ki = 0; ki < spec.strikes.size(); ++ki)
                    priceSums[ki] /= spec.rqmcReplicates;
            }

            for (std::size_t ki = 0; ki < spec.strikes.size(); ++ki)
                implied[ei][ki] = clampAndSolve(priceSums[ki], spec, ei, ki, params);
        }

        return VolSurface(spec.strikes, spec.expiries, implied);
    }

    static DriverCorrelationResult empiricalDriverCorrelation(
        const Parameters& params,
        double maturity,
        int steps,
        int paths,
        std::uint64_t seed = 42ULL)
    {
        if (paths <= 0)
            throw std::runtime_error("RoughBergomiModel: paths must be positive");

        Random rng(static_cast<unsigned int>(seed));
        std::vector<double> volIncrements;
        std::vector<double> spotIncrements;
        volIncrements.reserve(paths * steps);
        spotIncrements.reserve(paths * steps);

        for (int p = 0; p < paths; ++p) {
            PathResult path = simulatePath(params, 100.0, 0.0, maturity, steps, rng);
            for (int i = 0; i < steps; ++i) {
                volIncrements.push_back(path.volBrownianIncrements[i]);
                spotIncrements.push_back(path.spotBrownianIncrements[i]);
            }
        }

        DriverCorrelationResult result;
        result.targetRho = params.rho;
        result.empiricalRho = empiricalCorrelation(volIncrements, spotIncrements);
        result.samples = static_cast<int>(volIncrements.size());
        return result;
    }

private:
    static double clampAndSolve(
        double price,
        const SurfaceSpec& spec,
        std::size_t ei,
        std::size_t ki,
        const Parameters& params)
    {
        const double expiry = spec.expiries[ei];
        const double strike = spec.strikes[ki];
        const double intrinsic =
            std::max(spec.spot - strike * std::exp(-spec.rate * expiry), 0.0);
        const double upper = spec.spot;
        const double eps = 1e-8;

        price = std::max(price, intrinsic + eps);
        price = std::min(price, upper - eps);

        return ImpliedVolSolver::solve(
            price,
            spec.spot,
            strike,
            spec.rate,
            expiry,
            std::sqrt(params.xi0));
    }

    static double empiricalCorrelation(
        const std::vector<double>& x,
        const std::vector<double>& y)
    {
        if (x.size() != y.size() || x.empty())
            throw std::runtime_error("RoughBergomiModel: invalid correlation inputs");

        double meanX = std::accumulate(x.begin(), x.end(), 0.0) / x.size();
        double meanY = std::accumulate(y.begin(), y.end(), 0.0) / y.size();
        double cov = 0.0;
        double varX = 0.0;
        double varY = 0.0;
        for (std::size_t i = 0; i < x.size(); ++i) {
            const double dx = x[i] - meanX;
            const double dy = y[i] - meanY;
            cov += dx * dy;
            varX += dx * dx;
            varY += dy * dy;
        }
        if (varX <= 0.0 || varY <= 0.0)
            throw std::runtime_error("RoughBergomiModel: degenerate variance in correlation");
        return cov / std::sqrt(varX * varY);
    }

    static PricingResult summarizePayoffs(
        const std::vector<double>& payoffs,
        SamplingMethod method)
    {
        PricingResult result;
        result.price = std::accumulate(payoffs.begin(), payoffs.end(), 0.0) / payoffs.size();
        double sq = 0.0;
        for (std::size_t i = 0; i < payoffs.size(); ++i) {
            const double d = payoffs[i] - result.price;
            sq += d * d;
        }
        result.stdErr = std::sqrt(sq / std::max(1, static_cast<int>(payoffs.size()) - 1))
                      / std::sqrt(static_cast<double>(payoffs.size()));
        result.paths = static_cast<int>(payoffs.size());
        result.method = method;
        return result;
    }

    static std::uint64_t splitmix64(std::uint64_t x) {
        x += 0x9E3779B97F4A7C15ULL;
        x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
        x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
        return x ^ (x >> 31);
    }

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
        double q = 0.0;
        double r = 0.0;
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
