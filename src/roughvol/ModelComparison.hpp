#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <functional>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
#include "../calibration/ImpliedVolSolver.hpp"
#include "../payoffs/VanillaPayoffs.hpp"
#include "../pricing/BlackScholes.hpp"
#include "../pricing/EuropeanPricer.hpp"
#include "../integrators/EulerMaruyama.hpp"
#include "../simulators/MonteCarloSimulator.hpp"
#include "../simulators/MultiDimensionalMonteCarloSimulator.hpp"
#include "../stochasticProcess/Heston.hpp"
#include "../stochasticProcess/LocalVolProcess.hpp"
#include "../surface/LocalVolSurface.hpp"
#include "../surface/VolSurface.hpp"
#include "RoughBergomiCalibrator.hpp"
#include "RoughBergomiModel.hpp"

class ModelComparison {
public:
    struct HestonParameters {
        double v0;
        double kappa;
        double vbar;
        double xi;
        double rho;
    };

    struct SurfaceSpec {
        double spot;
        double rate;
        int steps;
        int paths;
        std::vector<double> strikes;
        std::vector<double> expiries;
    };

    struct BlackScholesFit {
        double sigma = 0.2;
        double rmse = 0.0;
        VolSurface surface = VolSurface({80.0, 100.0}, {0.5, 1.0}, {{0.2, 0.2}, {0.2, 0.2}});
    };

    struct HestonFit {
        HestonParameters params = {0.04, 1.0, 0.04, 0.4, -0.5};
        double rmse = 0.0;
        VolSurface surface = VolSurface({80.0, 100.0}, {0.5, 1.0}, {{0.2, 0.2}, {0.2, 0.2}});
    };

    static VolSurface generateBlackScholesSurface(
        double sigma,
        const SurfaceSpec& spec)
    {
        std::vector<std::vector<double>> vols(
            spec.expiries.size(),
            std::vector<double>(spec.strikes.size(), sigma));
        return VolSurface(spec.strikes, spec.expiries, vols);
    }

    static BlackScholesFit calibrateBlackScholes(
        const VolSurface& target,
        const SurfaceSpec& spec)
    {
        double bestSigma = 0.2;
        double bestRmse = std::numeric_limits<double>::infinity();
        for (double sigma = 0.05; sigma <= 1.0; sigma += 0.01) {
            const VolSurface surface = generateBlackScholesSurface(sigma, spec);
            const double error = RoughBergomiCalibrator::rmse(surface, target);
            if (error < bestRmse) {
                bestRmse = error;
                bestSigma = sigma;
            }
        }
        const VolSurface bestSurface = generateBlackScholesSurface(bestSigma, spec);
        return {bestSigma, bestRmse, bestSurface};
    }

    static VolSurface generateHestonSurface(
        const HestonParameters& params,
        const SurfaceSpec& spec)
    {
        std::vector<std::vector<double>> implied(
            spec.expiries.size(),
            std::vector<double>(spec.strikes.size(), 0.0));

        Heston heston(spec.rate, params.kappa, params.vbar, params.xi, params.rho);
        for (std::size_t ei = 0; ei < spec.expiries.size(); ++ei) {
            const double expiry = spec.expiries[ei];
            const double dt = expiry / std::max(1, spec.steps);
            const std::vector<double> x0 = {spec.spot, params.v0};
            const auto paths = MultiDimensionalMonteCarloSimulator::simulate(
                heston, x0, expiry, dt, spec.paths);

            std::vector<double> terminal(paths.size(), 0.0);
            for (std::size_t i = 0; i < paths.size(); ++i)
                terminal[i] = paths[i][0];

            for (std::size_t ki = 0; ki < spec.strikes.size(); ++ki) {
                CallPayoff payoff(spec.strikes[ki]);
                const auto priced = EuropeanPricer::priceWithError(terminal, payoff, spec.rate, expiry);
                const double intrinsic = std::max(spec.spot - spec.strikes[ki] * std::exp(-spec.rate * expiry), 0.0);
                const double upper = spec.spot;
                const double callPrice = std::max(intrinsic + 1e-8, std::min(priced.first, upper - 1e-8));
                implied[ei][ki] = ImpliedVolSolver::solve(
                    callPrice, spec.spot, spec.strikes[ki], spec.rate, expiry, std::sqrt(std::max(params.v0, 1e-4)));
            }
        }

        return VolSurface(spec.strikes, spec.expiries, implied);
    }

    static HestonFit calibrateHeston(
        const VolSurface& target,
        const SurfaceSpec& spec,
        const HestonParameters& initialGuess = {0.04, 1.0, 0.04, 0.5, -0.6},
        const HestonParameters& lowerBounds = {0.01, 0.2, 0.01, 0.1, -0.95},
        const HestonParameters& upperBounds = {0.20, 4.0, 0.20, 1.5, -0.05},
        int maxIterations = 18,
        double tolerance = 1e-4)
    {
        using Vec = std::array<double, 5>;
        auto toVec = [](const HestonParameters& p) {
            return Vec{{p.v0, p.kappa, p.vbar, p.xi, p.rho}};
        };
        auto clampValue = [](double x, double lo, double hi) {
            return std::max(lo, std::min(x, hi));
        };
        auto toParams = [&](const Vec& v) {
            HestonParameters p;
            p.v0 = clampValue(v[0], lowerBounds.v0, upperBounds.v0);
            p.kappa = clampValue(v[1], lowerBounds.kappa, upperBounds.kappa);
            p.vbar = clampValue(v[2], lowerBounds.vbar, upperBounds.vbar);
            p.xi = clampValue(v[3], lowerBounds.xi, upperBounds.xi);
            p.rho = clampValue(v[4], lowerBounds.rho, upperBounds.rho);
            return p;
        };

        const auto objective = [&](const Vec& v) {
            try {
                const VolSurface surface = generateHestonSurface(toParams(v), spec);
                return RoughBergomiCalibrator::rmse(surface, target);
            } catch (...) {
                return 1.0e6;
            }
        };

        Vec init = toVec(initialGuess);
        const int n = 5;
        std::vector<Vec> simplex(n + 1, init);
        for (int i = 0; i < n; ++i) {
            simplex[i + 1] = init;
            simplex[i + 1][i] += std::abs(init[i]) > 1e-6 ? 0.15 * std::abs(init[i]) : 0.02;
        }

        auto evaluateSimplex = [&](const std::vector<Vec>& points) {
            std::vector<std::pair<double, Vec>> scored;
            scored.reserve(points.size());
            for (const auto& point : points)
                scored.push_back({objective(point), point});
            std::sort(scored.begin(), scored.end(),
                [](const std::pair<double, Vec>& a, const std::pair<double, Vec>& b) { return a.first < b.first; });
            return scored;
        };

        for (int iter = 0; iter < maxIterations; ++iter) {
            auto scored = evaluateSimplex(simplex);
            for (int i = 0; i <= n; ++i)
                simplex[i] = scored[i].second;
            if (std::abs(scored.back().first - scored.front().first) < tolerance)
                break;

            Vec centroid = {};
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    centroid[j] += simplex[i][j] / n;

            Vec reflected = centroid;
            for (int j = 0; j < n; ++j)
                reflected[j] = centroid[j] + (centroid[j] - simplex.back()[j]);
            const double reflectedValue = objective(reflected);

            if (reflectedValue < scored.front().first) {
                Vec expanded = centroid;
                for (int j = 0; j < n; ++j)
                    expanded[j] = centroid[j] + 2.0 * (centroid[j] - simplex.back()[j]);
                simplex.back() = (objective(expanded) < reflectedValue) ? expanded : reflected;
            } else if (reflectedValue < scored[n - 1].first) {
                simplex.back() = reflected;
            } else {
                Vec contracted = centroid;
                for (int j = 0; j < n; ++j)
                    contracted[j] = centroid[j] + 0.5 * (simplex.back()[j] - centroid[j]);
                simplex.back() = (objective(contracted) < scored.back().first) ? contracted : simplex.front();
            }
        }

        const auto finalScored = evaluateSimplex(simplex);
        const HestonParameters bestParams = toParams(finalScored.front().second);
        const VolSurface bestSurface = generateHestonSurface(bestParams, spec);
        return {bestParams, finalScored.front().first, bestSurface};
    }

    static VolSurface generateLocalVolSurfaceFromImpliedSurface(
        const VolSurface& impliedSurface,
        const SurfaceSpec& spec)
    {
        LocalVolSurface localVol(impliedSurface, spec.spot, spec.rate);
        LocalVolProcess process(localVol, spec.rate);

        std::vector<std::vector<double>> implied(
            spec.expiries.size(),
            std::vector<double>(spec.strikes.size(), 0.0));

        const int localVolSteps = std::max(spec.steps * 4, 40);
        const int localVolPaths = std::max(spec.paths * 8, 2000);

        for (std::size_t ei = 0; ei < spec.expiries.size(); ++ei) {
            const double expiry = spec.expiries[ei];
            const double dt = expiry / std::max(1, localVolSteps);
            const auto terminal = MonteCarloSimulator<EulerMaruyama>::simulate(
                process, spec.spot, expiry, dt, localVolPaths);

            for (std::size_t ki = 0; ki < spec.strikes.size(); ++ki) {
                CallPayoff payoff(spec.strikes[ki]);
                const auto priced = EuropeanPricer::priceWithError(terminal, payoff, spec.rate, expiry);
                const double intrinsic = std::max(spec.spot - spec.strikes[ki] * std::exp(-spec.rate * expiry), 0.0);
                const double upper = spec.spot;
                const double callPrice = std::max(intrinsic + 1e-8, std::min(priced.first, upper - 1e-8));
                implied[ei][ki] = ImpliedVolSolver::solve(
                    callPrice,
                    spec.spot,
                    spec.strikes[ki],
                    spec.rate,
                    expiry,
                    localVol.smoothedImpliedVolAt(expiry, spec.strikes[ki]));
            }
        }

        return VolSurface(spec.strikes, spec.expiries, implied);
    }

    static void writeSurfaceSliceCsv(
        const std::string& path,
        const std::vector<std::string>& modelNames,
        const std::vector<VolSurface>& surfaces,
        const SurfaceSpec& spec,
        double smileExpiry,
        double termStrike)
    {
        std::ofstream out(path.c_str());
        out << "kind,model,x,y,implied_vol\n";
        for (std::size_t m = 0; m < surfaces.size(); ++m) {
            for (double K : spec.strikes)
                out << "smile," << modelNames[m] << "," << smileExpiry << "," << K << "," << surfaces[m].impliedVol(smileExpiry, K) << "\n";
            for (double T : spec.expiries)
                out << "term," << modelNames[m] << "," << termStrike << "," << T << "," << surfaces[m].impliedVol(T, termStrike) << "\n";
        }
    }

    static void writeRmseCsv(
        const std::string& path,
        const std::vector<std::string>& modelNames,
        const std::vector<double>& rmses)
    {
        std::ofstream out(path.c_str());
        out << "model,rmse\n";
        for (std::size_t i = 0; i < modelNames.size(); ++i)
            out << modelNames[i] << "," << rmses[i] << "\n";
    }
};
