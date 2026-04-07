#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include "../src/greeks/MultiAssetQMCGreeks.hpp"
#include "../src/payoffs/MultiAssetPayoff.hpp"
#include "../src/pricing/MultiAssetAnalytical.hpp"
#include "../src/pricing/MultiAssetPricer.hpp"
#include "../src/simulators/MultiAssetQMCSimulator.hpp"
#include "../src/stochasticProcess/CorrelatedGBM.hpp"

namespace {
double empiricalIncrementCorrelation(
    const std::vector<std::vector<double>>& terminalPrices,
    const std::vector<double>& spot0)
{
    double mean1 = 0.0;
    double mean2 = 0.0;
    for (const auto& path : terminalPrices) {
        mean1 += std::log(path[0] / spot0[0]);
        mean2 += std::log(path[1] / spot0[1]);
    }
    mean1 /= terminalPrices.size();
    mean2 /= terminalPrices.size();

    double cov = 0.0;
    double var1 = 0.0;
    double var2 = 0.0;
    for (const auto& path : terminalPrices) {
        const double x1 = std::log(path[0] / spot0[0]) - mean1;
        const double x2 = std::log(path[1] / spot0[1]) - mean2;
        cov += x1 * x2;
        var1 += x1 * x1;
        var2 += x2 * x2;
    }
    return cov / std::sqrt(var1 * var2);
}
}

int main() {
    const std::vector<double> mu = {0.05, 0.05};
    const std::vector<double> sigma = {0.20, 0.30};
    const double rho = 0.6;
    const std::vector<std::vector<double>> corr = {
        {1.0, rho},
        {rho, 1.0}
    };

    CorrelatedGBM process(mu, sigma, corr);
    const std::vector<double> S0 = {100.0, 100.0};
    const double r = 0.05;
    const double T = 1.0;
    const double dt = 0.01;
    const int paths = 4096;

    SpreadCallPayoff spread(0.0);
    const auto qmcPaths = MultiAssetQMCSimulator::simulate(
        process, S0, T, dt, paths, MultiAssetQMCSimulator::QMC,
        "docs/new-joe-kuo-6.21201.txt", true, false, 8, 42,
        MultiAssetQMCSimulator::Cholesky,
        MultiAssetQMCSimulator::AssetFirst);
    const auto qmcTimeFirstPaths = MultiAssetQMCSimulator::simulate(
        process, S0, T, dt, paths, MultiAssetQMCSimulator::QMC,
        "docs/new-joe-kuo-6.21201.txt", true, false, 8, 42,
        MultiAssetQMCSimulator::Cholesky,
        MultiAssetQMCSimulator::TimeFirst);
    const auto qmcHybridPaths = MultiAssetQMCSimulator::simulate(
        process, S0, T, dt, paths, MultiAssetQMCSimulator::QMC,
        "docs/new-joe-kuo-6.21201.txt", true, false, 8, 42,
        MultiAssetQMCSimulator::PcaBrownianBridgeHybrid);
    const auto qmcAdaptivePaths = MultiAssetQMCSimulator::simulate(
        process, S0, T, dt, paths, MultiAssetQMCSimulator::QMC,
        "docs/new-joe-kuo-6.21201.txt", true, false, 8, 42,
        MultiAssetQMCSimulator::Cholesky,
        MultiAssetQMCSimulator::AdaptiveVariance,
        &spread,
        16);
    const auto rqmcPaths = MultiAssetQMCSimulator::simulate(
        process, S0, T, dt, paths, MultiAssetQMCSimulator::RQMC,
        "docs/new-joe-kuo-6.21201.txt", true, true, 8);
    const auto rqmcHybridPaths = MultiAssetQMCSimulator::simulate(
        process, S0, T, dt, paths, MultiAssetQMCSimulator::RQMC,
        "docs/new-joe-kuo-6.21201.txt", true, true, 8, 42,
        MultiAssetQMCSimulator::PcaBrownianBridgeHybrid);

    const auto qmcPrice = MultiAssetPricer::price(qmcPaths.terminalPrices, spread, r, T);
    const auto qmcTimeFirstPrice = MultiAssetPricer::price(qmcTimeFirstPaths.terminalPrices, spread, r, T);
    const auto qmcHybridPrice = MultiAssetPricer::price(qmcHybridPaths.terminalPrices, spread, r, T);
    const auto qmcAdaptivePrice = MultiAssetPricer::price(qmcAdaptivePaths.terminalPrices, spread, r, T);
    const auto rqmcPrice = MultiAssetPricer::price(rqmcPaths.terminalPrices, spread, r, T);
    const auto rqmcHybridPrice = MultiAssetPricer::price(rqmcHybridPaths.terminalPrices, spread, r, T);
    const double analytical = MultiAssetAnalytical::margrabePrice(S0[0], S0[1], sigma[0], sigma[1], rho, T);

    const auto greeks = MultiAssetQMCGreeks::compute(
        process, S0, spread, r, T, dt, paths,
        MultiAssetQMCSimulator::QMC,
        "docs/new-joe-kuo-6.21201.txt",
        true, false, 8, 42,
        MultiAssetQMCSimulator::PcaBrownianBridgeHybrid,
        0.01, 0.01, 1e-3);

    const double analyticalDelta1 = (
        MultiAssetAnalytical::margrabePrice(S0[0] + 0.5, S0[1], sigma[0], sigma[1], rho, T) -
        MultiAssetAnalytical::margrabePrice(S0[0] - 0.5, S0[1], sigma[0], sigma[1], rho, T)) / 1.0;
    const double analyticalDelta2 = (
        MultiAssetAnalytical::margrabePrice(S0[0], S0[1] + 0.5, sigma[0], sigma[1], rho, T) -
        MultiAssetAnalytical::margrabePrice(S0[0], S0[1] - 0.5, sigma[0], sigma[1], rho, T)) / 1.0;

    const double empCorr = empiricalIncrementCorrelation(qmcPaths.terminalPrices, S0);

    std::cout << "=== Multi-Asset QMC Test ===\n\n";
    std::cout << "Sobol dimension used: " << qmcPaths.sobolDimension << "\n";
    std::cout << "QMC spread price:  " << std::fixed << std::setprecision(6) << qmcPrice.first << "\n";
    std::cout << "QMC time-first:    " << std::fixed << std::setprecision(6) << qmcTimeFirstPrice.first << "\n";
    std::cout << "QMC hybrid price:  " << std::fixed << std::setprecision(6) << qmcHybridPrice.first << "\n";
    std::cout << "QMC adaptive:      " << std::fixed << std::setprecision(6) << qmcAdaptivePrice.first << "\n";
    std::cout << "RQMC spread price: " << std::fixed << std::setprecision(6) << rqmcPrice.first << "\n";
    std::cout << "RQMC hybrid price: " << std::fixed << std::setprecision(6) << rqmcHybridPrice.first << "\n";
    std::cout << "Margrabe price:    " << std::fixed << std::setprecision(6) << analytical << "\n";
    std::cout << "Empirical corr:    " << std::fixed << std::setprecision(6) << empCorr
              << " (target " << rho << ")\n";
    std::cout << "QMC deltas:        [" << greeks.delta[0] << ", " << greeks.delta[1] << "]\n";
    std::cout << "Analytical approx: [" << analyticalDelta1 << ", " << analyticalDelta2 << "]\n";
    std::cout << "Adaptive top dimensions:";
    for (std::size_t i = 0; i < std::min<std::size_t>(6, qmcAdaptivePaths.orderedDimensions.size()); ++i)
        std::cout << " " << qmcAdaptivePaths.orderedDimensions[i];
    std::cout << "\n";

    assert(qmcPaths.sobolDimension == static_cast<int>(S0.size()) * static_cast<int>(T / dt));
    assert(std::abs(qmcPrice.first - analytical) < 0.35);
    assert(std::abs(qmcTimeFirstPrice.first - analytical) < 0.35);
    assert(std::abs(qmcHybridPrice.first - analytical) < 0.20);
    assert(std::abs(qmcAdaptivePrice.first - analytical) < 0.25);
    assert(std::abs(rqmcPrice.first - analytical) < 0.35);
    assert(std::abs(rqmcHybridPrice.first - analytical) < 0.20);
    assert(std::abs(empCorr - rho) < 0.08);
    assert(std::abs(greeks.delta[0] - analyticalDelta1) < 0.08);
    assert(std::abs(greeks.delta[1] - analyticalDelta2) < 0.08);
    assert(std::abs(qmcHybridPrice.first - analytical) <= std::abs(qmcPrice.first - analytical) + 0.05);
    assert(std::abs(qmcAdaptivePrice.first - analytical) <= std::abs(qmcTimeFirstPrice.first - analytical) + 0.10);
    assert(std::abs(rqmcHybridPrice.first - analytical) <= std::abs(rqmcPrice.first - analytical) + 0.05);

    std::cout << "\nMulti-asset Sobol/RQMC engine checks passed.\n";
    return 0;
}
