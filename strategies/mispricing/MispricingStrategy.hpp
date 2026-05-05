#pragma once
#include "../Strategy.hpp"
#include "../DeltaHedgedStraddleStrategy.hpp"
#include "MispricingScanner.hpp"
#include "../../surface/VolSurface.hpp"

// Runs a pair of delta-hedged positions based on relative value signals:
// long the undervalued option, short the overvalued one.
// Net PnL = short_pnl + long_pnl (convergence trade).
class MispricingStrategy : public Strategy {
public:
    struct Config {
        VolSurface market;
        VolSurface sviSurface;
        VolSurface rbergomiSurface;
        VolSurface localvolSurface;
        double spot;
        double hedgeInterval;
        double transactionCost;
        double minMispricing;   // threshold to trade
        int    maxLegs;
    };

    struct RunResult {
        BacktestResult aggregate;
        std::vector<MispricingScanner::RelativeValueTrade> trades;
        std::vector<MispricingScanner::MispricingPoint>    mispricings;
    };

    explicit MispricingStrategy(Config cfg) : config(std::move(cfg)) {}

    BacktestResult run(const std::vector<std::vector<double>>& paths,
                       double r, double T) override;

    RunResult runDetailed(const std::vector<std::vector<double>>& paths,
                          double r, double T);

    const char* name() const override {
        return "Mispricing RV Strategy";
    }

private:
    Config config;

    // Run a single-leg delta-hedged position and return per-path PnL
    std::vector<double> runLeg(
        const std::vector<std::vector<double>>& paths,
        double strike, double impliedVol, double r, double T,
        bool isShort) const;
};