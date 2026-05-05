#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>
#include "../surface/SVI.hpp"
#include "RoughBergomiModel.hpp"

class RoughBergomiCalibrator {
public:
    enum ObjectiveType {
        RawSurfaceRMSE,
        SviSmoothedSurfaceRMSE
    };

    struct CalibrationSpec {
        double spot;
        double rate;
        double maturity;
        int steps;
        int pathsPerEvaluation;
        int maxIterations;
        double tolerance;
        ObjectiveType objectiveType;
        RoughBergomiModel::Parameters initialGuess;
        RoughBergomiModel::Parameters lowerBounds;
        RoughBergomiModel::Parameters upperBounds;
        RoughBergomiModel::SamplingMethod objectiveSamplingMethod = RoughBergomiModel::MC;
        bool objectiveUseOwenScrambling = true;
        bool randomizeObjectiveSeed = false;
        std::uint64_t objectiveBaseSeed = 42ULL;
    };

    struct Result {
        RoughBergomiModel::Parameters parameters;
        double rmse;
        VolSurface fittedSurface;
    };

    static double rmse(
        const VolSurface& lhs,
        const VolSurface& rhs)
    {
        const std::vector<double>& strikes = lhs.getStrikes();
        const std::vector<double>& expiries = lhs.getExpiries();
        if (strikes.size() != rhs.getStrikes().size() || expiries.size() != rhs.getExpiries().size())
            throw std::runtime_error("RoughBergomiCalibrator: surfaces must share compatible grid sizes");

        double sq = 0.0;
        int count = 0;
        for (std::size_t i = 0; i < expiries.size(); ++i) {
            for (std::size_t j = 0; j < strikes.size(); ++j) {
                const double d = lhs.impliedVol(expiries[i], strikes[j])
                               - rhs.impliedVol(expiries[i], strikes[j]);
                sq += d * d;
                ++count;
            }
        }
        return std::sqrt(sq / std::max(1, count));
    }

    static Result calibrate(
        const VolSurface& target,
        const CalibrationSpec& spec)
    {
        validateSpec(spec);

        const SurfaceObjective objective = buildObjectiveSurface(target, spec);
        const RoughBergomiModel::SurfaceSpec surfaceSpec = buildSurfaceSpec(target, spec);
        std::uint64_t evaluationCounter = 0ULL;

        const ParameterVec optimum = nelderMead(
            [&](const ParameterVec& v) mutable {
                const RoughBergomiModel::Parameters p = toParameters(v, spec);
                try {
                    const int evalSeed = static_cast<int>(spec.objectiveBaseSeed + (spec.randomizeObjectiveSeed ? evaluationCounter++ : 0ULL));
                    const VolSurface modelSurface = RoughBergomiModel::generateVolSurface(
                        p,
                        surfaceSpec,
                        evalSeed,
                        spec.objectiveSamplingMethod,
                        spec.objectiveUseOwenScrambling);
                    return rmse(modelSurface, objective.surface);
                } catch (...) {
                    return 1.0e6;
                }
            },
            toVec(spec.initialGuess),
            spec.maxIterations,
            spec.tolerance);

        const RoughBergomiModel::Parameters bestParams = toParameters(optimum, spec);
        VolSurface bestSurface = objective.surface;
        double bestRmse = 1.0e6;
        try {
            bestSurface = RoughBergomiModel::generateVolSurface(
                bestParams,
                surfaceSpec,
                static_cast<int>(spec.objectiveBaseSeed),
                spec.objectiveSamplingMethod,
                spec.objectiveUseOwenScrambling);
            bestRmse = rmse(bestSurface, objective.surface);
        } catch (...) {
            bestSurface = objective.surface;
        }
        return {bestParams, bestRmse, bestSurface};
    }

private:
    typedef std::array<double, 4> ParameterVec;

    struct SurfaceObjective {
        VolSurface surface;
    };

    static void validateSpec(const CalibrationSpec& spec) {
        if (spec.pathsPerEvaluation <= 0 || spec.steps <= 0)
            throw std::runtime_error("RoughBergomiCalibrator: invalid evaluation grid");
        if (spec.maxIterations <= 0)
            throw std::runtime_error("RoughBergomiCalibrator: maxIterations must be positive");
        if (spec.tolerance <= 0.0)
            throw std::runtime_error("RoughBergomiCalibrator: tolerance must be positive");
    }

    static SurfaceObjective buildObjectiveSurface(
        const VolSurface& target,
        const CalibrationSpec& spec)
    {
        if (spec.objectiveType == RawSurfaceRMSE) {
            SurfaceObjective objective = {target};
            return objective;
        }

        const std::vector<double>& strikes = target.getStrikes();
        const std::vector<double>& expiries = target.getExpiries();
        std::vector<std::vector<double>> sviVols(
            expiries.size(),
            std::vector<double>(strikes.size(), 0.0));

        for (std::size_t ei = 0; ei < expiries.size(); ++ei) {
            const double T = expiries[ei];
            const double F = spec.spot * std::exp(spec.rate * T);
            std::vector<double> marketVols;
            marketVols.reserve(strikes.size());
            for (std::size_t ki = 0; ki < strikes.size(); ++ki)
                marketVols.push_back(target.impliedVol(T, strikes[ki]));

            const SVIParams fit = SVI::fit(strikes, marketVols, F, T);
            for (std::size_t ki = 0; ki < strikes.size(); ++ki)
            {
                try {
                    sviVols[ei][ki] = SVI::impliedVol(fit, strikes[ki], F, T);
                } catch (...) {
                    sviVols[ei][ki] = marketVols[ki];
                }
            }
        }

        SurfaceObjective objective = {VolSurface(strikes, expiries, sviVols)};
        return objective;
    }

    static RoughBergomiModel::SurfaceSpec buildSurfaceSpec(
        const VolSurface& target,
        const CalibrationSpec& spec)
    {
        RoughBergomiModel::SurfaceSpec surfaceSpec;
        surfaceSpec.spot = spec.spot;
        surfaceSpec.rate = spec.rate;
        surfaceSpec.maturity = spec.maturity;
        surfaceSpec.steps = spec.steps;
        surfaceSpec.paths = spec.pathsPerEvaluation;
        surfaceSpec.rqmcReplicates = 8;
        surfaceSpec.strikes = target.getStrikes();
        surfaceSpec.expiries = target.getExpiries();
        return surfaceSpec;
    }

    static ParameterVec toVec(const RoughBergomiModel::Parameters& p) {
        ParameterVec v;
        v[0] = p.xi0;
        v[1] = p.eta;
        v[2] = p.rho;
        v[3] = p.hurst;
        return v;
    }

    static double clampValue(double x, double lo, double hi) {
        return std::max(lo, std::min(x, hi));
    }

    static RoughBergomiModel::Parameters toParameters(
        const ParameterVec& v,
        const CalibrationSpec& spec)
    {
        RoughBergomiModel::Parameters p;
        p.xi0 = clampValue(v[0], spec.lowerBounds.xi0, spec.upperBounds.xi0);
        p.eta = clampValue(v[1], spec.lowerBounds.eta, spec.upperBounds.eta);
        p.rho = clampValue(v[2], spec.lowerBounds.rho, spec.upperBounds.rho);
        p.hurst = clampValue(v[3], spec.lowerBounds.hurst, spec.upperBounds.hurst);
        return p;
    }

    static ParameterVec nelderMead(
        const std::function<double(const ParameterVec&)>& f,
        const ParameterVec& init,
        int maxIter,
        double tol)
    {
        const int n = 4;
        std::vector<ParameterVec> simplex(n + 1, init);
        for (int i = 0; i < n; ++i) {
            simplex[i + 1] = init;
            simplex[i + 1][i] += std::abs(init[i]) > 1e-6 ? 0.1 * std::abs(init[i]) : 0.01;
        }

        auto evaluateSimplex = [&](const std::vector<ParameterVec>& points) {
            std::vector<std::pair<double, ParameterVec>> scored;
            scored.reserve(points.size());
            for (const auto& point : points)
                scored.push_back({f(point), point});
            std::sort(scored.begin(), scored.end(),
                [](const std::pair<double, ParameterVec>& a, const std::pair<double, ParameterVec>& b) {
                    return a.first < b.first;
                });
            return scored;
        };

        for (int iter = 0; iter < maxIter; ++iter) {
            std::vector<std::pair<double, ParameterVec>> scored = evaluateSimplex(simplex);
            for (int i = 0; i <= n; ++i)
                simplex[i] = scored[i].second;

            if (std::abs(scored.back().first - scored.front().first) < tol)
                break;

            ParameterVec centroid = {};
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    centroid[j] += simplex[i][j] / n;

            ParameterVec reflected = centroid;
            for (int j = 0; j < n; ++j)
                reflected[j] = centroid[j] + (centroid[j] - simplex.back()[j]);
            const double reflectedValue = f(reflected);

            if (reflectedValue < scored.front().first) {
                ParameterVec expanded = centroid;
                for (int j = 0; j < n; ++j)
                    expanded[j] = centroid[j] + 2.0 * (centroid[j] - simplex.back()[j]);
                const double expandedValue = f(expanded);
                simplex.back() = expandedValue < reflectedValue ? expanded : reflected;
            } else if (reflectedValue < scored[n - 1].first) {
                simplex.back() = reflected;
            } else {
                ParameterVec contracted = centroid;
                for (int j = 0; j < n; ++j)
                    contracted[j] = centroid[j] + 0.5 * (simplex.back()[j] - centroid[j]);
                const double contractedValue = f(contracted);
                if (contractedValue < scored.back().first) {
                    simplex.back() = contracted;
                } else {
                    for (int i = 1; i <= n; ++i)
                        for (int j = 0; j < n; ++j)
                            simplex[i][j] = simplex.front()[j]
                                          + 0.5 * (simplex[i][j] - simplex.front()[j]);
                }
            }
        }

        std::vector<std::pair<double, ParameterVec>> finalScored = evaluateSimplex(simplex);
        return finalScored.front().second;
    }
};
