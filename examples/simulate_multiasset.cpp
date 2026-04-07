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
double empiricalCorrelation(
    const std::vector<std::vector<double>>& terminalPrices,
    const std::vector<double>& spot0,
    int i,
    int j)
{
    double meanI = 0.0;
    double meanJ = 0.0;
    for (const auto& path : terminalPrices) {
        meanI += std::log(path[i] / spot0[i]);
        meanJ += std::log(path[j] / spot0[j]);
    }
    meanI /= terminalPrices.size();
    meanJ /= terminalPrices.size();

    double cov = 0.0;
    double varI = 0.0;
    double varJ = 0.0;
    for (const auto& path : terminalPrices) {
        const double xi = std::log(path[i] / spot0[i]) - meanI;
        const double xj = std::log(path[j] / spot0[j]) - meanJ;
        cov += xi * xj;
        varI += xi * xi;
        varJ += xj * xj;
    }
    return cov / std::sqrt(varI * varJ);
}

void printPriceLine(
    const std::string& label,
    const MultiAssetQMCSimulator::SimulationResult& simulation,
    const MultiAssetPayoff& payoff,
    double r,
    double T)
{
    const auto priced = MultiAssetPricer::price(simulation.terminalPrices, payoff, r, T);
    std::cout << "  " << std::left << std::setw(6) << label
              << " price=" << std::setw(10) << std::right << std::fixed << std::setprecision(4) << priced.first
              << " stdErr=" << std::setw(9) << priced.second
              << " sobolDim=" << simulation.sobolDimension
              << "\n";
}

double absError(double x, double ref) {
    return std::abs(x - ref);
}
}

int main() {
    const std::vector<double> mu = {0.05, 0.05, 0.05};
    const std::vector<double> sigma = {0.20, 0.25, 0.30};
    const std::vector<std::vector<double>> corrMatrix = {
        {1.00, 0.60, 0.25},
        {0.60, 1.00, 0.35},
        {0.25, 0.35, 1.00}
    };

    CorrelatedGBM cgbm(mu, sigma, corrMatrix);

    const std::vector<double> S0 = {100.0, 98.0, 105.0};
    const double r = 0.05;
    const double T = 1.0;
    const double dt = 0.02;
    const int paths = 4096;

    SpreadCallPayoff spread(0.0, 0, 1);
    BasketCallPayoff basket({0.4, 0.35, 0.25}, 100.0);

    std::cout << "Multi-Asset QMC Engine\n";
    std::cout << "  assets=" << S0.size()
              << " steps=" << static_cast<int>(T / dt)
              << " effective Sobol dimension=" << S0.size() * static_cast<int>(T / dt)
              << "\n";
    std::cout << "  correlation handled by Cholesky factor inside CorrelatedGBM\n\n";

    const auto mcSpread = MultiAssetQMCSimulator::simulate(
        cgbm, S0, T, dt, paths, MultiAssetQMCSimulator::MC);
    const auto qmcSpread = MultiAssetQMCSimulator::simulate(
        cgbm, S0, T, dt, paths, MultiAssetQMCSimulator::QMC,
        "docs/new-joe-kuo-6.21201.txt", true, false, 8, 42,
        MultiAssetQMCSimulator::Cholesky,
        MultiAssetQMCSimulator::AssetFirst);
    const auto qmcTimeFirstSpread = MultiAssetQMCSimulator::simulate(
        cgbm, S0, T, dt, paths, MultiAssetQMCSimulator::QMC,
        "docs/new-joe-kuo-6.21201.txt", true, false, 8, 42,
        MultiAssetQMCSimulator::Cholesky,
        MultiAssetQMCSimulator::TimeFirst);
    const auto qmcHybridSpread = MultiAssetQMCSimulator::simulate(
        cgbm, S0, T, dt, paths, MultiAssetQMCSimulator::QMC,
        "docs/new-joe-kuo-6.21201.txt", true, false, 8, 42,
        MultiAssetQMCSimulator::PcaBrownianBridgeHybrid);
    const auto qmcAdaptiveSpread = MultiAssetQMCSimulator::simulate(
        cgbm, S0, T, dt, paths, MultiAssetQMCSimulator::QMC,
        "docs/new-joe-kuo-6.21201.txt", true, false, 8, 42,
        MultiAssetQMCSimulator::Cholesky,
        MultiAssetQMCSimulator::AdaptiveVariance,
        &spread,
        16);
    const auto rqmcSpread = MultiAssetQMCSimulator::simulate(
        cgbm, S0, T, dt, paths, MultiAssetQMCSimulator::RQMC,
        "docs/new-joe-kuo-6.21201.txt", true, true, 8);
    const auto rqmcHybridSpread = MultiAssetQMCSimulator::simulate(
        cgbm, S0, T, dt, paths, MultiAssetQMCSimulator::RQMC,
        "docs/new-joe-kuo-6.21201.txt", true, true, 8, 42,
        MultiAssetQMCSimulator::PcaBrownianBridgeHybrid);

    const double margrabe = MultiAssetAnalytical::margrabePrice(
        S0[0], S0[1], sigma[0], sigma[1], corrMatrix[0][1], T);

    std::cout << "Spread Option (Margrabe benchmark)\n";
    printPriceLine("MC", mcSpread, spread, r, T);
    printPriceLine("QMC", qmcSpread, spread, r, T);
    printPriceLine("QMC+Time", qmcTimeFirstSpread, spread, r, T);
    printPriceLine("QMC+PCA", qmcHybridSpread, spread, r, T);
    printPriceLine("QMC+Adapt", qmcAdaptiveSpread, spread, r, T);
    printPriceLine("RQMC", rqmcSpread, spread, r, T);
    printPriceLine("RQMC+PCA", rqmcHybridSpread, spread, r, T);
    std::cout << "  Analytical benchmark=" << std::fixed << std::setprecision(4) << margrabe << "\n\n";

    std::cout << "Basket Call\n";
    printPriceLine("MC", mcSpread, basket, r, T);
    printPriceLine("QMC", qmcSpread, basket, r, T);
    printPriceLine("QMC+Time", qmcTimeFirstSpread, basket, r, T);
    printPriceLine("QMC+PCA", qmcHybridSpread, basket, r, T);
    printPriceLine("QMC+Adapt", qmcAdaptiveSpread, basket, r, T);
    printPriceLine("RQMC", rqmcSpread, basket, r, T);
    printPriceLine("RQMC+PCA", rqmcHybridSpread, basket, r, T);
    std::cout << "  Moment-match approx="
              << std::fixed << std::setprecision(4)
              << MultiAssetAnalytical::basketCallApprox(S0, {0.4, 0.35, 0.25}, sigma, corrMatrix, r, 100.0, T)
              << "\n\n";

    const auto plainSpreadPrice = MultiAssetPricer::price(qmcSpread.terminalPrices, spread, r, T).first;
    const auto timeFirstSpreadPrice = MultiAssetPricer::price(qmcTimeFirstSpread.terminalPrices, spread, r, T).first;
    const auto hybridSpreadPrice = MultiAssetPricer::price(qmcHybridSpread.terminalPrices, spread, r, T).first;
    const auto adaptiveSpreadPrice = MultiAssetPricer::price(qmcAdaptiveSpread.terminalPrices, spread, r, T).first;
    const auto plainBasketPrice = MultiAssetPricer::price(qmcSpread.terminalPrices, basket, r, T).first;
    const auto timeFirstBasketPrice = MultiAssetPricer::price(qmcTimeFirstSpread.terminalPrices, basket, r, T).first;
    const auto hybridBasketPrice = MultiAssetPricer::price(qmcHybridSpread.terminalPrices, basket, r, T).first;
    const auto adaptiveBasketPrice = MultiAssetPricer::price(qmcAdaptiveSpread.terminalPrices, basket, r, T).first;
    const double basketApprox = MultiAssetAnalytical::basketCallApprox(S0, {0.4, 0.35, 0.25}, sigma, corrMatrix, r, 100.0, T);

    std::cout << "Hybrid Construction Check\n";
    std::cout << "  spread |plain QMC - ref|=" << absError(plainSpreadPrice, margrabe)
              << " |time-first - ref|=" << absError(timeFirstSpreadPrice, margrabe)
              << " |hybrid QMC - ref|=" << absError(hybridSpreadPrice, margrabe) << "\n";
    std::cout << "  spread |adaptive QMC - ref|=" << absError(adaptiveSpreadPrice, margrabe) << "\n";
    std::cout << "  basket |plain QMC - approx|=" << absError(plainBasketPrice, basketApprox)
              << " |time-first - approx|=" << absError(timeFirstBasketPrice, basketApprox)
              << " |hybrid QMC - approx|=" << absError(hybridBasketPrice, basketApprox) << "\n\n";
    std::cout << "  basket |adaptive QMC - approx|=" << absError(adaptiveBasketPrice, basketApprox) << "\n\n";

    std::cout << "Adaptive ordering top coordinates";
    for (std::size_t i = 0; i < std::min<std::size_t>(8, qmcAdaptiveSpread.orderedDimensions.size()); ++i)
        std::cout << (i == 0 ? ": " : ", ") << qmcAdaptiveSpread.orderedDimensions[i];
    std::cout << "\n\n";

    std::cout << "Empirical Correlation Check (terminal log-returns)\n";
    std::cout << "  target rho(1,2)=" << corrMatrix[0][1]
              << " empirical=" << empiricalCorrelation(qmcSpread.terminalPrices, S0, 0, 1) << "\n";
    std::cout << "  target rho(1,3)=" << corrMatrix[0][2]
              << " empirical=" << empiricalCorrelation(qmcSpread.terminalPrices, S0, 0, 2) << "\n\n";

    const auto greeks = MultiAssetQMCGreeks::compute(
        cgbm, S0, basket, r, T, dt, paths,
        MultiAssetQMCSimulator::QMC,
        "docs/new-joe-kuo-6.21201.txt",
        true, false, 8, 42,
        MultiAssetQMCSimulator::PcaBrownianBridgeHybrid);

    std::cout << "Basket Greeks with QMC + PCA/Bridge Hybrid\n";
    std::cout << "  price=" << std::fixed << std::setprecision(4) << greeks.price << "\n";
    for (size_t i = 0; i < greeks.delta.size(); ++i) {
        std::cout << "  asset " << (i + 1)
                  << " delta=" << greeks.delta[i]
                  << " gamma=" << greeks.gamma[i]
                  << " vega=" << greeks.vega[i]
                  << "\n";
    }
    std::cout << "  rho=" << greeks.rho << "\n";

    return 0;
}
