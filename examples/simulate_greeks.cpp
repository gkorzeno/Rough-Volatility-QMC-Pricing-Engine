// examples/simulate_greeks.cpp
#include <iostream>
#include <iomanip>
#include "../src/greeks/Greeks.hpp"
#include "../src/greeks/AdjointGreeks.hpp"
#include "../src/greeks/FiniteDifferenceGreeks.hpp"
#include "../src/greeks/GreekEstimatorStudy.hpp"
#include "../src/greeks/PathwiseGreeks.hpp"
#include "../src/greeks/LikelihoodRatioGreeks.hpp"
#include "../src/payoffs/VanillaPayoffs.hpp"
#include "../src/pricing/BlackScholes.hpp"

int main() {
    double S0    = 100.0;
    double K     = 100.0;
    double r     = 0.05;
    double sigma = 0.2;
    double T     = 1.0;
    double dt    = 0.01;
    int    paths = 100000;

    CallPayoff call(K);

    auto bs = BlackScholes::greeks(S0, K, r, sigma, T);
    auto fd = FiniteDifferenceGreeks::compute(
                  call, S0, r, sigma, T, dt, paths);
    auto aad = AdjointGreeks::computeCall(
                  K, S0, r, sigma, T, dt, paths);
    auto pw = PathwiseGreeks::compute(
                  K, S0, r, sigma, T, dt, paths);
    auto lr = LikelihoodRatioGreeks::compute(
                  call, S0, r, sigma, T, dt, paths);

    std::cout << std::fixed << std::setprecision(4);
    std::cout << std::setw(12) << ""
              << std::setw(10) << "BS Exact"
              << std::setw(10) << "FD"
              << std::setw(10) << "AAD"
              << std::setw(10) << "Pathwise"
              << std::setw(10) << "LR\n";
    std::cout << std::string(62, '-') << "\n";

    auto row = [&](const std::string& name,
                   double b, double f, double a, double p, double l) {
        std::cout << std::setw(12) << name
                  << std::setw(10) << b
                  << std::setw(10) << f
                  << std::setw(10) << a
                  << std::setw(10) << p
                  << std::setw(10) << l << "\n";
    };

    row("Price",  bs.price, fd.price, aad.price, pw.price,  lr.price);
    row("Delta",  bs.delta, fd.delta, aad.delta, pw.delta,  lr.delta);
    row("Gamma",  bs.gamma, fd.gamma, 0.0,       0.0,       lr.gamma);
    row("Vega",   bs.vega,  fd.vega,  aad.vega,  pw.vega,   lr.vega);
    row("Theta",  bs.theta, fd.theta, 0.0,       0.0,       0.0);
    row("Rho",    bs.rho,   fd.rho,   aad.rho,   pw.rho,    0.0);

    // Digital call — LR shines here, pathwise breaks down
    std::cout << "\n-- Digital Call (K=100) --\n";
    DigitalCallPayoff digital(K);
    auto lr_dig = LikelihoodRatioGreeks::compute(
                      digital, S0, r, sigma, T, dt, paths);
    auto fd_dig = FiniteDifferenceGreeks::compute(
                      digital, S0, r, sigma, T, dt, paths);

    std::cout << std::setw(12) << ""
              << std::setw(10) << "FD"
              << std::setw(10) << "LR\n";
    std::cout << std::string(32, '-') << "\n";
    auto row2 = [&](const std::string& name, double f, double l) {
        std::cout << std::setw(12) << name
                  << std::setw(10) << f
                  << std::setw(10) << l << "\n";
    };
    row2("Price", fd_dig.price, lr_dig.price);
    row2("Delta", fd_dig.delta, lr_dig.delta);
    row2("Gamma", fd_dig.gamma, lr_dig.gamma);

    std::cout << "\n-- MC vs QMC vs RQMC estimator stability (Delta/Vega RMSE) --\n";
    for (const auto& sampler : {
            std::make_pair(GreekEstimatorStudy::MC, false),
            std::make_pair(GreekEstimatorStudy::QMC, false),
            std::make_pair(GreekEstimatorStudy::RQMC, true) }) {
        const auto summary = GreekEstimatorStudy::summarize(
            sampler.first, K, S0, r, sigma, T, dt, 2048, 6,
            "docs/new-joe-kuo-6.21201.txt", true, sampler.second, 8, 42);
        std::cout << GreekEstimatorStudy::samplerLabel(sampler.first, sampler.second) << "\n";
        std::cout << "  Pathwise         delta RMSE=" << summary.pathwise.deltaRmse
                  << " vega RMSE=" << summary.pathwise.vegaRmse << "\n";
        std::cout << "  LikelihoodRatio  delta RMSE=" << summary.likelihoodRatio.deltaRmse
                  << " vega RMSE=" << summary.likelihoodRatio.vegaRmse << "\n";
        std::cout << "  AAD              delta RMSE=" << summary.aad.deltaRmse
                  << " vega RMSE=" << summary.aad.vegaRmse << "\n";
    }

    return 0;
}
