#pragma once
#include "../Strategy.hpp"
#include "../DeltaHedgedStraddleStrategy.hpp"
#include "VolForecast.hpp"
#include "../../surface/VolSurface.hpp"
#include "../../roughvol/RoughBergomiModel.hpp"
#include <string>
#include <vector>

// Sells a delta-hedged straddle only when implied vol exceeds the rBergomi
// forecast by at least `entryThreshold` vol points.
// When the signal is absent the strategy stays flat (zero PnL for that path).
class VRPStrategy : public Strategy {
public:
    struct Config {
        RoughBergomiModel::Parameters modelParams;
        double impliedVol;         // market ATM IV
        double strike;
        double hedgeInterval;      // rebalancing frequency in years
        double entryThreshold;     // min VRP (IV - forecast) to enter, e.g. 0.02
        double transactionCost;
        int    forecastPaths;      // paths for vol forecast
        int    forecastSteps;
        std::string dirFile = "docs/new-joe-kuo-6.21201.txt";
    };

    struct Signal {
        double impliedVol;
        double forecastedVol;
        double vrp;
        bool   shouldTrade;
    };

    explicit VRPStrategy(Config cfg) : config(std::move(cfg)) {}

    // Compute the trading signal before running the backtest
    Signal computeSignal(double spot, double rate, double T) const {
        const auto fc = VolForecast::forecast(
            config.modelParams, spot, rate, T,
            config.forecastSteps, config.forecastPaths,
            RoughBergomiModel::RQMC, config.dirFile);

        const double vrp = VolForecast::vrp(config.impliedVol, fc.meanRealizedVol);
        return {config.impliedVol, fc.meanRealizedVol, vrp,
                vrp >= config.entryThreshold};
    }

    BacktestResult run(const std::vector<std::vector<double>>& paths,
                       double r, double T) override;

    const char* name() const override { return "VRP Strategy (rBergomi signal)"; }

    const Config& getConfig() const { return config; }

private:
    Config config;
};