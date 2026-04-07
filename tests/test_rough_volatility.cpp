#include <cassert>
#include <cmath>
#include <iostream>

#include "../src/core/Random.hpp"
#include "../src/roughvol/RoughVolatilityResearch.hpp"
#include "../src/surface/SVI.hpp"

int main() {
    Random rng(123);

    auto fBm = RoughVolatilityResearch::sampleFractionalBrownianMotion(1.0, 12, 0.15, rng);
    assert(static_cast<int>(fBm.size()) == 12);

    auto fracPath = RoughVolatilityResearch::simulateFractionalSDE(
        1.0,
        1.0,
        fBm,
        [](double x, double) { return -0.1 * x; },
        [](double, double) { return 0.25; });

    assert(fracPath.states.size() == 13);

    RoughVolatilityResearch::RoughBergomiParameters params = {0.04, 1.5, -0.7, 0.12};
    RoughVolatilityResearch::RoughSurfaceSpec spec;
    spec.spot = 100.0;
    spec.rate = 0.01;
    spec.maturity = 1.0;
    spec.steps = 12;
    spec.paths = 120;
    spec.rqmcReplicates = 8;
    spec.strikes = {90.0, 100.0, 110.0};
    spec.expiries = {0.5, 1.0};

    VolSurface surface = RoughVolatilityResearch::generateRoughBergomiSurface(params, spec, 123);
    const double atmHalf = surface.impliedVol(0.5, 100.0);
    const double atmOne = surface.impliedVol(1.0, 100.0);

    assert(atmHalf > 0.0);
    assert(atmOne > 0.0);

    auto smile = RoughVolatilityResearch::smileSlice(surface, 1.0, spec.strikes);
    double leftSkew = RoughVolatilityResearch::smileLeftSkewMetric(smile);
    assert(leftSkew > 0.0);

    auto term = RoughVolatilityResearch::termStructureSlice(surface, 100.0);
    assert(term.size() == spec.expiries.size());

    const double sviExpiry = 1.0;
    const double F = spec.spot * std::exp(spec.rate * sviExpiry);
    std::vector<double> marketVols;
    for (std::size_t i = 0; i < spec.strikes.size(); ++i)
        marketVols.push_back(surface.impliedVol(sviExpiry, spec.strikes[i]));
    SVIParams sviParams = SVI::fit(spec.strikes, marketVols, F, sviExpiry);
    double sviRmse = 0.0;
    for (std::size_t i = 0; i < spec.strikes.size(); ++i) {
        const double modelIv = SVI::impliedVol(sviParams, spec.strikes[i], F, sviExpiry);
        const double diff = modelIv - marketVols[i];
        sviRmse += diff * diff;
    }
    sviRmse = std::sqrt(sviRmse / spec.strikes.size());
    assert(sviRmse >= 0.0);

    RoughVolatilityResearch::CalibrationSpec calib;
    calib.spot = spec.spot;
    calib.rate = spec.rate;
    calib.maturity = spec.maturity;
    calib.steps = spec.steps;
    calib.pathsPerEvaluation = 20;
    calib.maxIterations = 4;
    calib.tolerance = 1e-4;
    calib.objectiveType = RoughBergomiCalibrator::SviSmoothedSurfaceRMSE;
    calib.initialGuess = {0.05, 1.1, -0.5, 0.18};
    calib.lowerBounds = {0.01, 0.20, -0.99, 0.02};
    calib.upperBounds = {0.20, 3.00, -0.05, 0.49};

    auto fit = RoughVolatilityResearch::calibrateRoughBergomi(surface, calib);
    assert(fit.rmse >= 0.0);
    assert(fit.parameters.hurst > 0.0 && fit.parameters.hurst < 0.5);

    RoughBergomiModel::PricingResult mcPrice = RoughBergomiModel::priceEuropeanCall(
        params, spec.spot, spec.rate, 1.0, spec.steps, 128, 100.0, RoughBergomiModel::MC);
    RoughBergomiModel::PricingResult qmcPrice = RoughBergomiModel::priceEuropeanCall(
        params, spec.spot, spec.rate, 1.0, spec.steps, 128, 100.0, RoughBergomiModel::QMC);
    RoughBergomiModel::PricingResult rqmcPrice = RoughBergomiModel::priceEuropeanCall(
        params, spec.spot, spec.rate, 1.0, spec.steps, 128, 100.0, RoughBergomiModel::RQMC,
        "docs/new-joe-kuo-6.21201.txt", true, 8);
    assert(mcPrice.price > 0.0);
    assert(qmcPrice.price > 0.0);
    assert(rqmcPrice.price > 0.0);
    assert(rqmcPrice.stdErr >= 0.0);

    RoughBergomiModel::DriverCorrelationResult corr =
        RoughBergomiModel::empiricalDriverCorrelation(params, 1.0, spec.steps, 300, 999ULL);
    assert(corr.samples == 300 * spec.steps);
    assert(std::abs(corr.empiricalRho - params.rho) < 0.10);

    auto variance = RoughVolatilityResearch::varianceAnalysis(
        params,
        spec,
        std::vector<int>{30, 60},
        180,
        3);
    assert(variance.size() == 2);
    assert(variance[0].rmse >= 0.0);
    assert(variance[1].rmse >= 0.0);
    assert(variance[0].rmseStdDev >= 0.0);
    assert(variance[1].rmseStdDev >= 0.0);

    std::cout << "rough-vol framework smoke test passed\n";
    std::cout << "ATM IV T=0.5: " << atmHalf << "\n";
    std::cout << "ATM IV T=1.0: " << atmOne << "\n";
    std::cout << "Left skew metric: " << leftSkew << "\n";
    std::cout << "SVI RMSE: " << sviRmse << "\n";
    std::cout << "Calibration RMSE: " << fit.rmse << "\n";
    std::cout << "MC/QMC/RQMC prices: "
              << mcPrice.price << " / "
              << qmcPrice.price << " / "
              << rqmcPrice.price << "\n";
    std::cout << "Target vs empirical driver rho: "
              << params.rho << " / "
              << corr.empiricalRho << "\n";
    return 0;
}
