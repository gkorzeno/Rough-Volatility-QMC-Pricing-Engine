#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include "RoughPathStorage.hpp"
#include "LSMCBasis.hpp"
#include "LSMCRegressor.hpp"
#include "../../payoffs/Payoff.hpp"

// Longstaff-Schwartz Monte Carlo for American options under rough Bergomi.
//
// The non-Markovian state (S, v, Y) is used as the regression state at each
// exercise date. This correctly captures path-dependence: two paths at the
// same spot price but different variance histories can have very different
// continuation values under rough vol.
//
// Algorithm:
//   1. Simulate full rBergomi paths storing (S, v, Y) at every step.
//   2. Initialise cashflow = payoff(S_T) at maturity.
//   3. Backward: at each step j (from T-1 to 1):
//      a. Find ITM paths (intrinsic > 0).
//      b. Build basis features from (S_j, v_j, Y_j) for ITM paths.
//      c. Regress discounted future cashflow onto basis → beta_j.
//      d. Compare intrinsic vs predicted continuation; update cashflow.
//   4. Price = mean of discounted cashflows.
class AmericanPricer {
public:
    struct Config {
        int    degree          = 3;
        LSMCBasis::Type basis  = LSMCBasis::Laguerre;
        bool   includeVol      = true;   // use v in regression state
        bool   includeVolterra = true;   // use Y in regression state
        double regularisation  = 1e-6;
        int    minITM          = 50;     // min ITM paths to attempt regression
    };

    struct Result {
        double price;
        double stdErr;
        double earlyExercisePremium;  // vs European
        double europeanPrice;
        int    nPaths;
        int    nSteps;
        // Regression R^2 by exercise date (diagnostic)
        std::vector<double> rSquaredByStep;
        // Exercise boundary: mean S at exercise at each step
        std::vector<double> exerciseBoundary;
    };

    static Result price(
        const RoughBergomiModel::Parameters& modelParams,
        const Payoff&  payoff,
        double spot,
        double rate,
        double maturity,
        int    nSteps,
        int    nPaths,
        const Config& cfg,
        unsigned int seed = 42)
    {
        // ── 1. Simulate full paths ────────────────────────────────────────
        auto data = RoughPathStorage::simulate(
            modelParams, spot, rate, maturity, nSteps, nPaths, seed);

        const double dt   = maturity / nSteps;
        const double disc = std::exp(-rate * dt);  // per-step discount factor

        // ── 2. Initialise cashflow at maturity ────────────────────────────
        std::vector<double> cashflow(nPaths);
        for (int i = 0; i < nPaths; ++i)
            cashflow[i] = payoff(data.S[i][nSteps]);

        // ── 3. Diagnostics storage ────────────────────────────────────────
        std::vector<double> rSquaredByStep(nSteps, 0.0);
        std::vector<double> exerciseBoundary(nSteps, 0.0);

        // ── 4. Backward induction ─────────────────────────────────────────
        for (int step = nSteps - 1; step >= 1; --step) {
            const double tau = data.times[step];

            // Cross-section at this step
            std::vector<double> S_j, v_j, Y_j;
            RoughPathStorage::crossSection(data, step, S_j, v_j, Y_j);

            // Identify ITM paths
            std::vector<int> itmIdx;
            itmIdx.reserve(nPaths / 2);
            for (int i = 0; i < nPaths; ++i)
                if (payoff(S_j[i]) > 0.0)
                    itmIdx.push_back(i);

            if (static_cast<int>(itmIdx.size()) < cfg.minITM) {
                // Not enough ITM paths — discount and continue
                for (int i = 0; i < nPaths; ++i)
                    cashflow[i] *= disc;
                continue;
            }

            // ── 4a. Build regression inputs for ITM paths ─────────────────
            std::vector<std::vector<double>> X_itm;
            std::vector<double>              y_itm;
            X_itm.reserve(itmIdx.size());
            y_itm.reserve(itmIdx.size());

            for (int i : itmIdx) {
                X_itm.push_back(LSMCBasis::build(
                    S_j[i], v_j[i], Y_j[i],
                    spot, modelParams.xi0,
                    cfg.degree, cfg.basis,
                    cfg.includeVol, cfg.includeVolterra));
                // Response: discounted future cashflow
                y_itm.push_back(disc * cashflow[i]);
            }

            // ── 4b. Regress ───────────────────────────────────────────────
            auto fit = LSMCRegressor::fit(X_itm, y_itm, cfg.regularisation);
            rSquaredByStep[step] = fit.rSquared;

            // ── 4c. Exercise decision for ITM paths ───────────────────────
            double sumExerciseS = 0.0;
            int    nExercised   = 0;

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
                    // Exercise: replace discounted future cashflow with intrinsic
                    cashflow[i] = intrinsic;
                    sumExerciseS += S_j[i];
                    ++nExercised;
                } else {
                    // Hold: discount the future cashflow
                    cashflow[i] *= disc;
                }
            }

            // Discount OTM paths (they cannot exercise here)
            for (int i = 0; i < nPaths; ++i) {
                if (payoff(S_j[i]) <= 0.0)
                    cashflow[i] *= disc;
            }

            exerciseBoundary[step] = (nExercised > 0)
                ? sumExerciseS / nExercised : 0.0;
        }

        // ── 5. Discount back to t=0 ───────────────────────────────────────
        // cashflow[i] is already at t=1 after the backward loop above;
        // discount one more step to t=0.
        for (int i = 0; i < nPaths; ++i)
            cashflow[i] *= disc;

        // ── 6. Price and std error ────────────────────────────────────────
        const double mean = std::accumulate(
            cashflow.begin(), cashflow.end(), 0.0) / nPaths;

        double sq = 0.0;
        for (double c : cashflow) sq += (c - mean) * (c - mean);
        const double se = std::sqrt(sq / std::max(1, nPaths - 1)) /
                          std::sqrt(static_cast<double>(nPaths));

        // European price for comparison (no early exercise)
        double euroMean = 0.0;
        const double totalDisc = std::exp(-rate * maturity);
        for (int i = 0; i < nPaths; ++i)
            euroMean += payoff(data.S[i][nSteps]);
        euroMean = totalDisc * euroMean / nPaths;

        return {mean, se, mean - euroMean, euroMean,
                nPaths, nSteps, rSquaredByStep, exerciseBoundary};
    }

    static Result price(
        const RoughBergomiModel::Parameters& modelParams,
        const Payoff&  payoff,
        double spot,
        double rate,
        double maturity,
        int    nSteps,
        int    nPaths,
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
            Config{},
            seed);
    }
};
