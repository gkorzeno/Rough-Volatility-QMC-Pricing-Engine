#pragma once
#include "../Strategy.hpp"
#include "../DeltaHedgedStraddleStrategy.hpp"
#include "SkewAnalyzer.hpp"
#include "../../surface/VolSurface.hpp"
#include <vector>

// Implements calendar spread and risk reversal trades derived from
// SVI surface analysis and rBergomi term structure dynamics.
class TermStructureStrategy : public Strategy {
public:
    struct Config {
        VolSurface surface;
        double spot;
        double rate;
        double hedgeInterval;
        double transactionCost;
        double skewThreshold;
        double termSpreadThreshold;
        int    maxSignals;
    };

    struct RunResult {
        BacktestResult aggregate;
        std::vector<SkewAnalyzer::Signal> signals;
        std::vector<SkewAnalyzer::SkewMetrics> skewByExpiry;
    };

    explicit TermStructureStrategy(Config cfg) : config(std::move(cfg)) {}

    BacktestResult run(const std::vector<std::vector<double>>& paths,
                       double r, double T) override;

    RunResult runDetailed(const std::vector<std::vector<double>>& paths,
                          double r, double T);

    const char* name() const override { return "Term Structure / Skew Strategy"; }

private:
    Config config;
};