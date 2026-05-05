#pragma once
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#include "Greeks.hpp"
#include "../payoffs/MultiAssetPayoff.hpp"
#include "../pricing/MultiAssetPricer.hpp"
#include "../simulators/MultiAssetQMCSimulator.hpp"
#include "../stochasticProcess/CorrelatedGBM.hpp"

struct MultiAssetGreeks {
    double price = 0.0;
    std::vector<double> delta;
    std::vector<double> gamma;
    std::vector<double> vega;
    double rho = 0.0;
};

class MultiAssetQMCGreeks {
public:
    static MultiAssetGreeks compute(
        const CorrelatedGBM& process,
        const std::vector<double>& S0,
        const MultiAssetPayoff& payoff,
        double r,
        double T,
        double dt,
        int paths,
        MultiAssetQMCSimulator::SamplingMethod method = MultiAssetQMCSimulator::QMC,
        const std::string& directionFile = "docs/new-joe-kuo-6.21201.txt",
        bool useBrownianBridge = true,
        bool useOwenScrambling = false,
        int rqmcReplicates = 8,
        std::uint64_t seed = 42,
        MultiAssetQMCSimulator::PathConstruction construction = MultiAssetQMCSimulator::Cholesky,
        double relativeSpotShift = 0.01,
        double absoluteVolShift = 0.01,
        double rateShift = 1e-3)
    {
        MultiAssetGreeks greeks;
        const int assets = static_cast<int>(S0.size());
        greeks.delta.assign(assets, 0.0);
        greeks.gamma.assign(assets, 0.0);
        greeks.vega.assign(assets, 0.0);

        const auto baseTerminal = simulate(
            process, S0, T, dt, paths, method, directionFile,
            useBrownianBridge, useOwenScrambling, rqmcReplicates, seed, construction);
        greeks.price = MultiAssetPricer::price(baseTerminal.terminalPrices, payoff, r, T).first;

        for (int i = 0; i < assets; ++i) {
            const double dS = std::max(1e-4, std::abs(S0[i]) * relativeSpotShift);

            auto spotUp = S0;
            auto spotDown = S0;
            spotUp[i] += dS;
            spotDown[i] = std::max(1e-8, spotDown[i] - dS);

            const auto upTerminal = simulate(
                process, spotUp, T, dt, paths, method, directionFile,
                useBrownianBridge, useOwenScrambling, rqmcReplicates, seed, construction);
            const auto downTerminal = simulate(
                process, spotDown, T, dt, paths, method, directionFile,
                useBrownianBridge, useOwenScrambling, rqmcReplicates, seed, construction);

            const double upPrice = MultiAssetPricer::price(upTerminal.terminalPrices, payoff, r, T).first;
            const double downPrice = MultiAssetPricer::price(downTerminal.terminalPrices, payoff, r, T).first;

            greeks.delta[i] = (upPrice - downPrice) / (spotUp[i] - spotDown[i]);
            greeks.gamma[i] = (upPrice - 2.0 * greeks.price + downPrice)
                / std::pow(0.5 * (spotUp[i] - spotDown[i]), 2.0);
        }

        const auto baseMu = process.drifts();
        const auto baseSigma = process.volatilities();
        const auto baseCorr = process.correlation();

        for (int i = 0; i < assets; ++i) {
            auto sigmaUp = baseSigma;
            auto sigmaDown = baseSigma;
            sigmaUp[i] += absoluteVolShift;
            sigmaDown[i] = std::max(1e-6, sigmaDown[i] - absoluteVolShift);

            CorrelatedGBM upProcess(baseMu, sigmaUp, baseCorr);
            CorrelatedGBM downProcess(baseMu, sigmaDown, baseCorr);

            const auto upTerminal = simulate(
                upProcess, S0, T, dt, paths, method, directionFile,
                useBrownianBridge, useOwenScrambling, rqmcReplicates, seed, construction);
            const auto downTerminal = simulate(
                downProcess, S0, T, dt, paths, method, directionFile,
                useBrownianBridge, useOwenScrambling, rqmcReplicates, seed, construction);

            const double upPrice = MultiAssetPricer::price(upTerminal.terminalPrices, payoff, r, T).first;
            const double downPrice = MultiAssetPricer::price(downTerminal.terminalPrices, payoff, r, T).first;
            greeks.vega[i] = (upPrice - downPrice) / (sigmaUp[i] - sigmaDown[i]);
        }

        std::vector<double> muUp(assets, r + rateShift);
        std::vector<double> muDown(assets, r - rateShift);
        CorrelatedGBM rhoUpProcess(muUp, baseSigma, baseCorr);
        CorrelatedGBM rhoDownProcess(muDown, baseSigma, baseCorr);

        const auto rhoUpTerminal = simulate(
            rhoUpProcess, S0, T, dt, paths, method, directionFile,
            useBrownianBridge, useOwenScrambling, rqmcReplicates, seed, construction);
        const auto rhoDownTerminal = simulate(
            rhoDownProcess, S0, T, dt, paths, method, directionFile,
            useBrownianBridge, useOwenScrambling, rqmcReplicates, seed, construction);

        const double rhoUpPrice = MultiAssetPricer::price(rhoUpTerminal.terminalPrices, payoff, r + rateShift, T).first;
        const double rhoDownPrice = MultiAssetPricer::price(rhoDownTerminal.terminalPrices, payoff, r - rateShift, T).first;
        greeks.rho = (rhoUpPrice - rhoDownPrice) / (2.0 * rateShift);

        return greeks;
    }

private:
    static MultiAssetQMCSimulator::SimulationResult simulate(
        const CorrelatedGBM& process,
        const std::vector<double>& S0,
        double T,
        double dt,
        int paths,
        MultiAssetQMCSimulator::SamplingMethod method,
        const std::string& directionFile,
        bool useBrownianBridge,
        bool useOwenScrambling,
        int rqmcReplicates,
        std::uint64_t seed,
        MultiAssetQMCSimulator::PathConstruction construction)
    {
        return MultiAssetQMCSimulator::simulate(
            process, S0, T, dt, paths, method, directionFile,
            useBrownianBridge, useOwenScrambling, rqmcReplicates, seed, construction);
    }
};
