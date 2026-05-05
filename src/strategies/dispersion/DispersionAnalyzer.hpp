#pragma once
#include <vector>
#include <cmath>
#include <string>
#include <stdexcept>
#include "../../pricing/MultiAssetPricer.hpp"
#include "../../pricing/BlackScholes.hpp"
#include "../../simulators/MultiAssetQMCSimulator.hpp"
#include "../../stochasticProcess/CorrelatedGBM.hpp"
#include "../../payoffs/MultiAssetPayoff.hpp"

class DispersionAnalyzer {
public:
    struct DispersionMetrics {
        double indexImpliedVol;       // IV of the basket option
        double avgConstituentIV;      // avg IV of individual options
        double impliedCorrelation;    // backed out from index vol
        double modelCorrelation;      // from CorrelatedGBM
        double dispersionSpread;      // avgConstituentIV - indexImpliedVol
        bool   dispersionOpportunity; // true when spread is large
    };

    // Compute implied correlation from index and constituent IVs.
    // Uses the standard identity:
    //   sigma_index^2 = sum_i sum_j w_i w_j sigma_i sigma_j rho_ij
    // For equal weights and uniform pairwise correlation rho:
    //   sigma_I^2 = sigma_avg^2 * [rho + (1-rho)/n]
    // => rho = (sigma_I^2/sigma_avg^2 - 1/n) / (1 - 1/n)
    // Replace impliedCorrelation with the correct general formula:
    static double impliedCorrelation(
        double indexIV,
        const std::vector<double>& constituentIVs,
        const std::vector<double>& weights)
    {
        const int n = static_cast<int>(constituentIVs.size());
        if (n < 2) throw std::runtime_error("Need at least 2 constituents");

        // Normalise weights
        double sumW = 0.0;
        for (double w : weights) sumW += w;

        // Weighted average constituent variance and cross terms
        // sigma_I^2 = sum_i w_i^2 sigma_i^2  +  rho * sum_{i!=j} w_i w_j sigma_i sigma_j
        // Solve for rho:
        //   rho = (sigma_I^2 - sum_i (w_i/W)^2 sigma_i^2)
        //         / (sum_{i!=j} (w_i/W)(w_j/W) sigma_i sigma_j)

        double diagSum  = 0.0;  // sum w_i^2 sigma_i^2
        double crossSum = 0.0;  // sum_{i!=j} w_i w_j sigma_i sigma_j
        for (int i = 0; i < n; ++i) {
            const double wi = weights[i] / sumW;
            diagSum += wi * wi * constituentIVs[i] * constituentIVs[i];
            for (int j = 0; j < n; ++j) {
                if (i == j) continue;
                const double wj = weights[j] / sumW;
                crossSum += wi * wj * constituentIVs[i] * constituentIVs[j];
            }
        }

        if (crossSum < 1e-12) return 0.0;

        const double indexVar = indexIV * indexIV;
        const double rho = (indexVar - diagSum) / crossSum;
        return std::max(-1.0, std::min(1.0, rho));
    }
    // Compare model implied correlation with market implied correlation
    static DispersionMetrics analyze(
        double indexIV,
        const std::vector<double>& constituentIVs,
        const std::vector<double>& weights,
        double modelCorrelation,
        double dispersionThreshold = 0.02)
    {
        double avgIV = 0.0, sumW = 0.0;
        for (std::size_t i = 0; i < constituentIVs.size(); ++i) {
            avgIV += weights[i] * constituentIVs[i];
            sumW  += weights[i];
        }
        avgIV /= sumW;

        const double implCorr = impliedCorrelation(indexIV, constituentIVs, weights);
        const double spread   = avgIV - indexIV;

        return {indexIV, avgIV, implCorr, modelCorrelation,
                spread, spread >= dispersionThreshold};
    }
};