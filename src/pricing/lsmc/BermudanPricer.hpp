#pragma once
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include "RoughPathStorage.hpp"
#include "LSMCBasis.hpp"
#include "LSMCRegressor.hpp"
#include "../../payoffs/Payoff.hpp"

// Bermudan option: exercise only allowed on a fixed set of dates.
// Under rough vol this is distinct from American in a meaningful way:
// fewer exercise opportunities means the rBergomi term structure of vol
// affects which dates are most valuable for early exercise.
class BermudanPricer {
public:
    using Config = AmericanPricer::Config;

    struct Result {
        double price;
        double stdErr;
        double earlyExercisePremium;
        double europeanPrice;
        // Per exercise-date diagnostics
        struct ExerciseDateResult {
            double time;
            int    step;
            double rSquared;
            double exerciseBoundary;  // mean S at exercise
            double exerciseFraction;  // fraction of ITM paths exercising
        };
        std::vector<ExerciseDateResult> byDate;
    };

    // exerciseTimes: sorted list of allowed exercise times in (0, maturity]
    // Must align with the simulation time grid.
    static Result price(
        const RoughBergomiModel::Parameters& modelParams,
        const Payoff&  payoff,
        double spot,
        double rate,
        double maturity,
        int    nSteps,
        int    nPaths,
        const std::vector<double>& exerciseTimes,
        const Config& cfg,
        unsigned int seed = 42)
    {
        if (exerciseTimes.empty())
            throw std::runtime_error("BermudanPricer: exerciseTimes must be non-empty");
        if (exerciseTimes.back() > maturity + 1e-10)
            throw std::runtime_error("BermudanPricer: exercise times must be <= maturity");

        // Snap exercise times to nearest grid steps
        const double dt = maturity / nSteps;
        std::vector<int> exerciseSteps;
        for (double t : exerciseTimes) {
            int step = static_cast<int>(std::round(t / dt));
            step = std::max(1, std::min(step, nSteps));
            exerciseSteps.push_back(step);
        }
        // Deduplicate and sort descending (for backward induction)
        std::sort(exerciseSteps.begin(), exerciseSteps.end());
        exerciseSteps.erase(
            std::unique(exerciseSteps.begin(), exerciseSteps.end()),
            exerciseSteps.end());
        std::vector<int> stepsDesc(exerciseSteps.rbegin(), exerciseSteps.rend());

        // ── Simulate ──────────────────────────────────────────────────────
        auto data = RoughPathStorage::simulate(
            modelParams, spot, rate, maturity, nSteps, nPaths, seed);

        const double discStep = std::exp(-rate * dt);

        // Initialise cashflow at maturity
        std::vector<double> cashflow(nPaths);
        for (int i = 0; i < nPaths; ++i)
            cashflow[i] = payoff(data.S[i][nSteps]);

        // Track which step each path last touched (for diagnostics)
        std::vector<Result::ExerciseDateResult> byDate;

        // Build a set for O(1) lookup of exercise steps
        std::vector<bool> isExerciseStep(nSteps + 1, false);
        for (int s : exerciseSteps) isExerciseStep[s] = true;

        // ── Backward induction ────────────────────────────────────────────
        for (int step = nSteps - 1; step >= 1; --step) {
            // Discount all cashflows one step
            for (int i = 0; i < nPaths; ++i)
                cashflow[i] *= discStep;

            if (!isExerciseStep[step]) continue;

            std::vector<double> S_j, v_j, Y_j;
            RoughPathStorage::crossSection(data, step, S_j, v_j, Y_j);

            std::vector<int> itmIdx;
            for (int i = 0; i < nPaths; ++i)
                if (payoff(S_j[i]) > 0.0)
                    itmIdx.push_back(i);

            Result::ExerciseDateResult edr;
            edr.time = data.times[step];
            edr.step = step;
            edr.rSquared = 0.0;
            edr.exerciseBoundary = 0.0;
            edr.exerciseFraction = 0.0;

            if (static_cast<int>(itmIdx.size()) < cfg.minITM) {
                byDate.push_back(edr);
                continue;
            }

            // Build regression
            std::vector<std::vector<double>> X_itm;
            std::vector<double>              y_itm;
            for (int i : itmIdx) {
                X_itm.push_back(LSMCBasis::build(
                    S_j[i], v_j[i], Y_j[i],
                    spot, modelParams.xi0,
                    cfg.degree, cfg.basis,
                    cfg.includeVol, cfg.includeVolterra));
                y_itm.push_back(cashflow[i]);
            }

            auto fit = LSMCRegressor::fit(X_itm, y_itm, cfg.regularisation);
            edr.rSquared = fit.rSquared;

            double sumExS = 0.0;
            int    nEx    = 0;

            for (int i : itmIdx) {
                const double intrinsic = payoff(S_j[i]);
                const auto   features  = LSMCBasis::build(
                    S_j[i], v_j[i], Y_j[i],
                    spot, modelParams.xi0,
                    cfg.degree, cfg.basis,
                    cfg.includeVol, cfg.includeVolterra);
                const double continuation = LSMCRegressor::predict(
                    fit.coefficients, features);

                if (intrinsic >= continuation) {
                    cashflow[i] = intrinsic;
                    sumExS += S_j[i];
                    ++nEx;
                }
            }

            edr.exerciseBoundary = (nEx > 0) ? sumExS / nEx : 0.0;
            edr.exerciseFraction = (itmIdx.empty()) ? 0.0
                : static_cast<double>(nEx) / itmIdx.size();
            byDate.push_back(edr);
        }

        // Sort by date ascending for output
        std::reverse(byDate.begin(), byDate.end());

        // ── Price ─────────────────────────────────────────────────────────
        const double mean = std::accumulate(
            cashflow.begin(), cashflow.end(), 0.0) / nPaths;
        double sq = 0.0;
        for (double c : cashflow) sq += (c - mean) * (c - mean);
        const double se = std::sqrt(sq / std::max(1, nPaths - 1)) /
                          std::sqrt(static_cast<double>(nPaths));

        const double totalDisc = std::exp(-rate * maturity);
        double euroMean = 0.0;
        for (int i = 0; i < nPaths; ++i)
            euroMean += payoff(data.S[i][nSteps]);
        euroMean = totalDisc * euroMean / nPaths;

        return {mean, se, mean - euroMean, euroMean, byDate};
    }

    static Result price(
        const RoughBergomiModel::Parameters& modelParams,
        const Payoff&  payoff,
        double spot,
        double rate,
        double maturity,
        int    nSteps,
        int    nPaths,
        const std::vector<double>& exerciseTimes,
        unsigned int seed = 42)
    {
        return price(
            modelParams,
            payoff,
            spot,
            rate,
            maturity,
            nSteps,
            nPaths,
            exerciseTimes,
            Config{},
            seed);
    }
};

// Make BermudanPricer visible to AmericanPricer.hpp by completing the forward ref
#include "AmericanPricer.hpp"
