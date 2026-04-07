#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include "../payoffs/MultiAssetPayoff.hpp"
#include "../core/BrownianBridge.hpp"
#include "../core/Random.hpp"
#include "../core/SobolSequence.hpp"
#include "../integrators/MultiDimensionalEulerMaruyama.hpp"
#include "../stochasticProcess/CorrelatedGBM.hpp"
#include "../stochasticProcess/MultiDimensionalProcess.hpp"

class MultiAssetQMCSimulator {
public:
    enum SamplingMethod {
        MC,
        QMC,
        RQMC
    };

    enum PathConstruction {
        Cholesky,
        PcaBrownianBridgeHybrid
    };

    enum DimensionOrdering {
        AssetFirst,
        TimeFirst,
        AdaptiveVariance
    };

    struct SimulationResult {
        std::vector<std::vector<double>> terminalPrices;
        int assets = 0;
        int steps = 0;
        int sobolDimension = 0;
        int replicates = 1;
        bool usedBrownianBridge = true;
        bool usedOwenScrambling = false;
        PathConstruction construction = Cholesky;
        DimensionOrdering ordering = AssetFirst;
        std::vector<int> orderedDimensions;
        std::vector<double> dimensionImportance;
    };

    static SimulationResult simulate(
        const MultiDimensionalProcess& process,
        const std::vector<double>& x0,
        double T,
        double dt,
        int paths,
        SamplingMethod method = QMC,
        const std::string& directionFile = "docs/new-joe-kuo-6.21201.txt",
        bool useBrownianBridge = true,
        bool useOwenScrambling = false,
        int rqmcReplicates = 8,
        std::uint64_t seed = 42,
        PathConstruction construction = Cholesky,
        DimensionOrdering ordering = AssetFirst,
        const MultiAssetPayoff* adaptivePayoff = nullptr,
        int adaptivePilotPaths = 32)
    {
        if (x0.empty())
            throw std::runtime_error("MultiAssetQMCSimulator: x0 must be non-empty");
        if (T <= 0.0 || dt <= 0.0)
            throw std::runtime_error("MultiAssetQMCSimulator: T and dt must be positive");

        const int steps = static_cast<int>(T / dt + 1e-12);
        if (steps <= 0)
            throw std::runtime_error("MultiAssetQMCSimulator: invalid time grid");

        const int assets = static_cast<int>(x0.size());
        const int sobolDimension = assets * steps;

        SimulationResult result;
        result.assets = assets;
        result.steps = steps;
        result.sobolDimension = sobolDimension;
        result.usedBrownianBridge = useBrownianBridge;
        result.usedOwenScrambling = useOwenScrambling;
        result.replicates = (method == RQMC) ? rqmcReplicates : 1;
        result.construction = construction;
        result.ordering = ordering;
        result.terminalPrices.resize(paths);

        if (construction == PcaBrownianBridgeHybrid && !useBrownianBridge)
            throw std::runtime_error("MultiAssetQMCSimulator: PCA hybrid requires Brownian bridge ON");

        PcaData pcaData;
        if (construction == PcaBrownianBridgeHybrid) {
            const auto* correlated = dynamic_cast<const CorrelatedGBM*>(&process);
            if (correlated == nullptr)
                throw std::runtime_error("MultiAssetQMCSimulator: PCA hybrid currently requires CorrelatedGBM");
            pcaData = buildPcaData(*correlated);
        }

        if (ordering == AdaptiveVariance && adaptivePayoff == nullptr)
            throw std::runtime_error("MultiAssetQMCSimulator: adaptive ordering requires a payoff");

        if (ordering == AdaptiveVariance && adaptivePilotPaths <= 0)
            throw std::runtime_error("MultiAssetQMCSimulator: adaptivePilotPaths must be positive");

        std::vector<double> dimensionImportance;
        std::vector<int> orderedDimensions = buildDimensionOrder(
            process, x0, T, dt, assets, steps, useBrownianBridge, construction,
            ordering, adaptivePayoff, adaptivePilotPaths, seed, pcaData, &dimensionImportance);
        result.orderedDimensions = orderedDimensions;
        result.dimensionImportance = dimensionImportance;

        if (method == MC) {
            Random rng(static_cast<unsigned int>(seed));
            for (int path = 0; path < paths; ++path) {
                if (construction == PcaBrownianBridgeHybrid) {
                    auto correlatedNormals = mcPcaHybridNormals(assets, steps, rng, dt, pcaData);
                    result.terminalPrices[path] = evolveCorrelatedGbmPath(process, x0, T, dt, correlatedNormals);
                } else {
                    auto normals = mcNormals(assets, steps, rng, useBrownianBridge, dt);
                    result.terminalPrices[path] = evolvePath(process, x0, T, dt, normals);
                }
            }
            return result;
        }

        if (method == QMC) {
            SobolSequence sobol(sobolDimension, directionFile);
            for (int path = 0; path < paths; ++path) {
                auto uniforms = sobol.next();
                if (construction == PcaBrownianBridgeHybrid) {
                    auto correlatedNormals = sobolPcaHybridNormals(
                        assets, steps, uniforms, dt, pcaData, orderedDimensions);
                    result.terminalPrices[path] = evolveCorrelatedGbmPath(process, x0, T, dt, correlatedNormals);
                } else {
                    auto normals = sobolNormals(
                        assets, steps, uniforms, useBrownianBridge, dt, orderedDimensions);
                    result.terminalPrices[path] = evolvePath(process, x0, T, dt, normals);
                }
            }
            return result;
        }

        if (rqmcReplicates <= 1 || paths % rqmcReplicates != 0)
            throw std::runtime_error("MultiAssetQMCSimulator: RQMC paths must be divisible by replicates and replicates>1");

        const int pathsPerReplicate = paths / rqmcReplicates;
        for (int rep = 0; rep < rqmcReplicates; ++rep) {
            const std::uint64_t repSeed = seed + 104729ULL * static_cast<std::uint64_t>(rep + 1);
            SobolSequence sobol = useOwenScrambling
                ? SobolSequence(sobolDimension, directionFile, repSeed, 12)
                : SobolSequence(sobolDimension, directionFile);

            for (int i = 0; i < pathsPerReplicate; ++i) {
                auto uniforms = sobol.next();
                if (!useOwenScrambling)
                    applyDigitalShift(uniforms, repSeed);

                const int pathIndex = rep * pathsPerReplicate + i;
                if (construction == PcaBrownianBridgeHybrid) {
                    auto correlatedNormals = sobolPcaHybridNormals(
                        assets, steps, uniforms, dt, pcaData, orderedDimensions);
                    result.terminalPrices[pathIndex] = evolveCorrelatedGbmPath(process, x0, T, dt, correlatedNormals);
                } else {
                    auto normals = sobolNormals(
                        assets, steps, uniforms, useBrownianBridge, dt, orderedDimensions);
                    result.terminalPrices[pathIndex] = evolvePath(process, x0, T, dt, normals);
                }
            }
        }

        return result;
    }

private:
    struct PcaData {
        std::vector<double> eigenvalues;
        std::vector<std::vector<double>> eigenvectors;
        std::vector<int> order;
        std::vector<std::vector<double>> loadings;
    };

    static double probit(double u)
    {
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
        double q, r;
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

    static std::vector<std::vector<double>> mcNormals(
        int assets,
        int steps,
        Random& rng,
        bool useBrownianBridge,
        double dt)
    {
        std::vector<std::vector<double>> normals(steps, std::vector<double>(assets, 0.0));
        for (int asset = 0; asset < assets; ++asset) {
            std::vector<double> pathNormals(steps);
            for (int step = 0; step < steps; ++step)
                pathNormals[step] = rng.normal();

            assignTemporalPath(normals, asset, pathNormals, useBrownianBridge, dt);
        }
        return normals;
    }

    static std::vector<std::vector<double>> sobolNormals(
        int assets,
        int steps,
        const std::vector<double>& uniforms,
        bool useBrownianBridge,
        double dt,
        const std::vector<int>& orderedDimensions)
    {
        if (static_cast<int>(uniforms.size()) != assets * steps)
            throw std::runtime_error("MultiAssetQMCSimulator: uniform vector size mismatch");

        std::vector<std::vector<double>> rawPaths(assets, std::vector<double>(steps, 0.0));
        for (int q = 0; q < assets * steps; ++q) {
            const int flat = orderedDimensions.empty() ? q : orderedDimensions[q];
            const int asset = flat / steps;
            const int step = flat % steps;
            rawPaths[asset][step] = probit(uniforms[q]);
        }

        std::vector<std::vector<double>> normals(steps, std::vector<double>(assets, 0.0));
        for (int asset = 0; asset < assets; ++asset) {
            assignTemporalPath(normals, asset, rawPaths[asset], useBrownianBridge, dt);
        }
        return normals;
    }

    static std::vector<std::vector<double>> mcPcaHybridNormals(
        int assets,
        int steps,
        Random& rng,
        double dt,
        const PcaData& pcaData)
    {
        std::vector<std::vector<double>> factorNormals(steps, std::vector<double>(assets, 0.0));
        for (int rank = 0; rank < assets; ++rank) {
            std::vector<double> pathNormals(steps);
            for (int step = 0; step < steps; ++step)
                pathNormals[step] = rng.normal();

            auto dW = BrownianBridge::build(pathNormals, dt);
            const double invSqrtDt = 1.0 / std::sqrt(dt);
            for (int step = 0; step < steps; ++step)
                factorNormals[step][rank] = dW[step] * invSqrtDt;
        }
        return combinePcaFactors(factorNormals, pcaData);
    }

    static std::vector<std::vector<double>> sobolPcaHybridNormals(
        int assets,
        int steps,
        const std::vector<double>& uniforms,
        double dt,
        const PcaData& pcaData,
        const std::vector<int>& orderedDimensions)
    {
        if (static_cast<int>(uniforms.size()) != assets * steps)
            throw std::runtime_error("MultiAssetQMCSimulator: uniform vector size mismatch");

        std::vector<std::vector<double>> rawFactorPaths(assets, std::vector<double>(steps, 0.0));
        for (int q = 0; q < assets * steps; ++q) {
            const int flat = orderedDimensions.empty() ? q : orderedDimensions[q];
            const int rank = flat / steps;
            const int step = flat % steps;
            rawFactorPaths[rank][step] = probit(uniforms[q]);
        }

        std::vector<std::vector<double>> factorNormals(steps, std::vector<double>(assets, 0.0));
        for (int rank = 0; rank < assets; ++rank) {
            const int factorIndex = pcaData.order[rank];
            auto dW = BrownianBridge::build(rawFactorPaths[rank], dt);
            const double invSqrtDt = 1.0 / std::sqrt(dt);
            for (int step = 0; step < steps; ++step)
                factorNormals[step][factorIndex] = dW[step] * invSqrtDt;
        }
        return combinePcaFactors(factorNormals, pcaData);
    }

    static void assignTemporalPath(
        std::vector<std::vector<double>>& normals,
        int asset,
        const std::vector<double>& pathNormals,
        bool useBrownianBridge,
        double dt)
    {
        if (useBrownianBridge) {
            auto dW = BrownianBridge::build(pathNormals, dt);
            const double invSqrtDt = 1.0 / std::sqrt(dt);
            for (size_t step = 0; step < pathNormals.size(); ++step)
                normals[step][asset] = dW[step] * invSqrtDt;
        } else {
            for (size_t step = 0; step < pathNormals.size(); ++step)
                normals[step][asset] = pathNormals[step];
        }
    }

    static std::vector<std::vector<double>> combinePcaFactors(
        const std::vector<std::vector<double>>& factorNormals,
        const PcaData& pcaData)
    {
        const int steps = static_cast<int>(factorNormals.size());
        const int assets = static_cast<int>(factorNormals.front().size());
        std::vector<std::vector<double>> correlated(steps, std::vector<double>(assets, 0.0));
        for (int step = 0; step < steps; ++step) {
            for (int asset = 0; asset < assets; ++asset) {
                double value = 0.0;
                for (int factor = 0; factor < assets; ++factor)
                    value += pcaData.loadings[asset][factor] * factorNormals[step][factor];
                correlated[step][asset] = value;
            }
        }
        return correlated;
    }

    static std::vector<double> evolvePath(
        const MultiDimensionalProcess& process,
        const std::vector<double>& x0,
        double T,
        double dt,
        const std::vector<std::vector<double>>& normals)
    {
        const int steps = static_cast<int>(T / dt + 1e-12);
        std::vector<double> x = x0;
        double t = 0.0;
        for (int step = 0; step < steps; ++step) {
            x = MultiDimensionalEulerMaruyama::stepWithZ(process, x, t, dt, normals[step]);
            t += dt;
        }
        return x;
    }

    static std::vector<double> evolveCorrelatedGbmPath(
        const MultiDimensionalProcess& process,
        const std::vector<double>& x0,
        double T,
        double dt,
        const std::vector<std::vector<double>>& correlatedNormals)
    {
        const auto* cgbm = dynamic_cast<const CorrelatedGBM*>(&process);
        if (cgbm == nullptr)
            throw std::runtime_error("MultiAssetQMCSimulator: correlated GBM path construction requires CorrelatedGBM");

        const int steps = static_cast<int>(T / dt + 1e-12);
        const auto& mu = cgbm->drifts();
        const auto& sigma = cgbm->volatilities();
        const int assets = static_cast<int>(x0.size());
        const double sqrtDt = std::sqrt(dt);

        std::vector<double> x = x0;
        for (int step = 0; step < steps; ++step) {
            for (int asset = 0; asset < assets; ++asset)
                x[asset] += mu[asset] * x[asset] * dt
                    + sigma[asset] * x[asset] * correlatedNormals[step][asset] * sqrtDt;
        }
        return x;
    }

    static void applyDigitalShift(std::vector<double>& uniforms, std::uint64_t seed)
    {
        for (size_t i = 0; i < uniforms.size(); ++i) {
            const std::uint64_t mixed = seed
                ^ (0x9E3779B97F4A7C15ULL + static_cast<std::uint64_t>(i) * 0xBF58476D1CE4E5B9ULL);
            const double shift = static_cast<double>((mixed >> 11) & ((1ULL << 53) - 1))
                / static_cast<double>(1ULL << 53);
            double shifted = uniforms[i] + shift;
            shifted -= std::floor(shifted);
            if (shifted <= 0.0) shifted = 1e-15;
            if (shifted >= 1.0) shifted = 1.0 - 1e-15;
            uniforms[i] = shifted;
        }
    }

    static PcaData buildPcaData(const CorrelatedGBM& process)
    {
        const auto corr = process.correlation();
        const auto eigen = jacobiEigenDecomposition(corr);
        const int n = static_cast<int>(corr.size());

        std::vector<int> order(n);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            return eigen.first[a] > eigen.first[b];
        });

        PcaData data;
        data.eigenvalues = eigen.first;
        data.eigenvectors = eigen.second;
        data.order = order;
        data.loadings.assign(n, std::vector<double>(n, 0.0));

        for (int rank = 0; rank < n; ++rank) {
            const int factor = order[rank];
            const double scale = std::sqrt(std::max(0.0, eigen.first[factor]));
            for (int asset = 0; asset < n; ++asset)
                data.loadings[asset][factor] = eigen.second[asset][factor] * scale;
        }
        return data;
    }

    static std::vector<int> buildDimensionOrder(
        const MultiDimensionalProcess& process,
        const std::vector<double>& x0,
        double T,
        double dt,
        int assets,
        int steps,
        bool useBrownianBridge,
        PathConstruction construction,
        DimensionOrdering ordering,
        const MultiAssetPayoff* adaptivePayoff,
        int adaptivePilotPaths,
        std::uint64_t seed,
        const PcaData& pcaData,
        std::vector<double>* importanceOut)
    {
        const int dim = assets * steps;
        if (importanceOut != nullptr)
            importanceOut->assign(dim, 0.0);

        if (ordering == AssetFirst)
            return assetFirstOrder(assets, steps);
        if (ordering == TimeFirst)
            return timeFirstOrder(assets, steps);

        std::vector<double> importance(dim, 0.0);
        Random rng(static_cast<unsigned int>(seed + 7919));
        std::vector<double> uniforms(dim, 0.0);
        std::vector<double> modified(dim, 0.0);

        for (int pilot = 0; pilot < adaptivePilotPaths; ++pilot) {
            for (int i = 0; i < dim; ++i) {
                double u = rng.uniform();
                if (u <= 0.0) u = 1e-15;
                if (u >= 1.0) u = 1.0 - 1e-15;
                uniforms[i] = u;
            }

            const std::vector<double> baseTerminal = (construction == PcaBrownianBridgeHybrid)
                ? evolveCorrelatedGbmPath(
                    process, x0, T, dt,
                    sobolPcaHybridNormals(assets, steps, uniforms, dt, pcaData, assetFirstOrder(assets, steps)))
                : evolvePath(
                    process, x0, T, dt,
                    sobolNormals(assets, steps, uniforms, useBrownianBridge, dt, assetFirstOrder(assets, steps)));
            const double basePayoff = (*adaptivePayoff)(baseTerminal);

            for (int j = 0; j < dim; ++j) {
                modified = uniforms;
                modified[j] = 0.5;
                const std::vector<double> altTerminal = (construction == PcaBrownianBridgeHybrid)
                    ? evolveCorrelatedGbmPath(
                        process, x0, T, dt,
                        sobolPcaHybridNormals(assets, steps, modified, dt, pcaData, assetFirstOrder(assets, steps)))
                    : evolvePath(
                        process, x0, T, dt,
                        sobolNormals(assets, steps, modified, useBrownianBridge, dt, assetFirstOrder(assets, steps)));
                const double diff = basePayoff - (*adaptivePayoff)(altTerminal);
                importance[j] += diff * diff;
            }
        }

        std::vector<int> order = assetFirstOrder(assets, steps);
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            return importance[a] > importance[b];
        });
        if (importanceOut != nullptr)
            *importanceOut = importance;
        return order;
    }

    static std::vector<int> assetFirstOrder(int assets, int steps)
    {
        std::vector<int> order(assets * steps);
        std::iota(order.begin(), order.end(), 0);
        return order;
    }

    static std::vector<int> timeFirstOrder(int assets, int steps)
    {
        std::vector<int> order;
        order.reserve(assets * steps);
        for (int step = 0; step < steps; ++step)
            for (int asset = 0; asset < assets; ++asset)
                order.push_back(asset * steps + step);
        return order;
    }

    static std::pair<std::vector<double>, std::vector<std::vector<double>>> jacobiEigenDecomposition(
        const std::vector<std::vector<double>>& input)
    {
        const int n = static_cast<int>(input.size());
        std::vector<std::vector<double>> a = input;
        std::vector<std::vector<double>> v(n, std::vector<double>(n, 0.0));
        for (int i = 0; i < n; ++i)
            v[i][i] = 1.0;

        for (int iter = 0; iter < 100 * n * n; ++iter) {
            int p = 0;
            int q = 1;
            double maxOffDiag = 0.0;
            for (int i = 0; i < n; ++i) {
                for (int j = i + 1; j < n; ++j) {
                    if (std::abs(a[i][j]) > maxOffDiag) {
                        maxOffDiag = std::abs(a[i][j]);
                        p = i;
                        q = j;
                    }
                }
            }

            if (maxOffDiag < 1e-12)
                break;

            const double app = a[p][p];
            const double aqq = a[q][q];
            const double apq = a[p][q];
            const double tau = (aqq - app) / (2.0 * apq);
            const double t = ((tau >= 0.0) ? 1.0 : -1.0)
                / (std::abs(tau) + std::sqrt(1.0 + tau * tau));
            const double c = 1.0 / std::sqrt(1.0 + t * t);
            const double s = t * c;

            for (int k = 0; k < n; ++k) {
                if (k == p || k == q)
                    continue;
                const double aik = a[k][p];
                const double akq = a[k][q];
                a[k][p] = c * aik - s * akq;
                a[p][k] = a[k][p];
                a[k][q] = s * aik + c * akq;
                a[q][k] = a[k][q];
            }

            a[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
            a[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
            a[p][q] = 0.0;
            a[q][p] = 0.0;

            for (int k = 0; k < n; ++k) {
                const double vip = v[k][p];
                const double viq = v[k][q];
                v[k][p] = c * vip - s * viq;
                v[k][q] = s * vip + c * viq;
            }
        }

        std::vector<double> eigenvalues(n);
        for (int i = 0; i < n; ++i)
            eigenvalues[i] = a[i][i];
        return {eigenvalues, v};
    }
};
