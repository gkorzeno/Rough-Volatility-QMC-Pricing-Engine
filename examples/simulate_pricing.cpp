// examples/simulate_pricing.cpp
#include <iostream>
#include "../src/stochasticProcess/GeometricBrownianMotion.hpp"
#include "../src/simulators/MonteCarloSimulator.hpp"
#include "../src/simulators/AntitheticSimulator.hpp"
#include "../src/simulators/QMCSimulator.hpp"
#include "../src/integrators/EulerMaruyama.hpp"
#include "../src/integrators/Milstein.hpp"
#include "../src/payoffs/VanillaPayoffs.hpp"
#include "../src/pricing/EuropeanPricer.hpp"
#include "../src/pricing/BlackScholes.hpp"
#include "../src/pricing/ControlVariatePricer.hpp"

int main() {
    double S0    = 100.0;
    double K     = 100.0;
    double r     = 0.05;
    double sigma = 0.2;
    double T     = 1.0;
    double dt    = 0.01;
    int    paths = 10000;

    GeometricBrownianMotion gbm(r, sigma);
    CallPayoff call(K);

    double bsPrice = BlackScholes::callPrice(S0, K, r, sigma, T);

    // Plain MC
    auto res_em = MonteCarloSimulator<EulerMaruyama>::simulate(gbm, S0, T, dt, paths);
    auto [price_em, err_em] = EuropeanPricer::priceWithError(res_em, call, r, T);

    // Milstein
    auto res_mil = MonteCarloSimulator<Milstein>::simulate(gbm, S0, T, dt, paths);
    auto [price_mil, err_mil] = EuropeanPricer::priceWithError(res_mil, call, r, T);

    // Antithetic
    auto res_anti = AntitheticSimulator<EulerMaruyama>::simulate(gbm, S0, T, dt, paths);
    auto [price_anti, err_anti] = EuropeanPricer::priceWithError(res_anti, call, r, T);

    // Control variate
    auto [price_cv, ratio] = ControlVariatePricer::priceWithDiagnostics(
        res_em, call, r, T, S0);

    // QMC
    // auto res_qmc = QMCSimulator<EulerMaruyama>::simulate(gbm, S0, T, dt, paths);
    auto res_qmc = QMCSimulator<EulerMaruyama>::simulate(
        gbm, S0, T, dt, paths, "docs/new-joe-kuo-6.21201.txt");
    auto [price_qmc, err_qmc] = EuropeanPricer::priceWithError(res_qmc, call, r, T);

    std::cout << "Black-Scholes analytical: " << bsPrice          << "\n";
    std::cout << "Euler-Maruyama MC:        " << price_em
              << " ± " << err_em                                   << "\n";
    std::cout << "Milstein MC:              " << price_mil
              << " ± " << err_mil                                  << "\n";
    std::cout << "Antithetic MC:            " << price_anti
              << " ± " << err_anti                                 << "\n";
    std::cout << "Control Variate:          " << price_cv
              << " (var reduction: " << (1.0 - ratio)*100 << "%)" << "\n";
    std::cout << "Quasi-MC (Sobol):         " << price_qmc
              << " ± " << err_qmc                                  << "\n";

    return 0;
}