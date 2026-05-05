#pragma once
#include "../Strategy.hpp"
#include "../DeltaHedgedStraddleStrategy.hpp"
#include "DispersionAnalyzer.hpp"
#include "../../stochasticProcess/CorrelatedGBM.hpp"
#include "../../simulators/MultiAssetQMCSimulator.hpp"
#include "../../pricing/MultiAssetPricer.hpp"
#include "../../payoffs/MultiAssetPayoff.hpp"
#include <vector>
#include <numeric>
#include <algorithm>

class DispersionStrategy : public Strategy {
public:
    struct Config {
        double               indexStrike;
        double               indexImpliedVol;
        std::vector<double>  basketWeights;      // must sum to 1
        std::vector<double>  constituentStrikes;
        std::vector<double>  constituentIVs;
        CorrelatedGBM        process;
        std::vector<double>  spotVec;
        double hedgeInterval;
        double transactionCost;
        double entryThreshold;  // min avgConstituentIV - indexIV to enter
        std::string dirFile = "docs/new-joe-kuo-6.21201.txt";
    };

    struct RunResult {
        BacktestResult aggregate;
        DispersionAnalyzer::DispersionMetrics metrics;
        bool traded;
        double realizedIndexVol;       // ex-post realised vol of basket
        double realizedConstituentVol; // mean ex-post realised vol of constituents
    };

    explicit DispersionStrategy(Config cfg) : config(std::move(cfg)) {}

    // Strategy interface: generates its own correlated paths internally
    BacktestResult run(const std::vector<std::vector<double>>& /*ignored*/,
                       double r, double T) override {
        return runDispersion(r, T, 10000).aggregate;
    }

    // Full run with configurable path count and diagnostics
    RunResult runDispersion(double r, double T, int nPaths,
                            std::uint64_t seed = 42);

    const char* name() const override { return "Dispersion Trading Strategy"; }

private:
    Config config;

    // Compute realised vol from a set of 1D paths (annualised)
    static double realisedVol(
        const std::vector<std::vector<double>>& paths, double T)
    {
        if (paths.empty() || paths[0].size() < 2) return 0.0;
        const int nSteps = static_cast<int>(paths[0].size()) - 1;
        const double dt = T / nSteps;
        double sumVar = 0.0;
        for (const auto& path : paths) {
            double pathVar = 0.0;
            for (int j = 1; j <= nSteps; ++j) {
                const double logRet = std::log(path[j] / path[j-1]);
                pathVar += logRet * logRet;
            }
            sumVar += pathVar / (nSteps * dt);
        }
        return std::sqrt(sumVar / paths.size());
    }
};