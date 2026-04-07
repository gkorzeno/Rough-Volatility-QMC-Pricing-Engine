#include <iostream>

#include "../src/stochasticProcess/GeometricBrownianMotion.hpp"
#include "../src/simulators/MonteCarloSimulator.hpp"
#include "../src/integrators/EulerMaruyama.hpp"
#include "../src/pricing/EuropeanPricer.hpp"
#include "../src/pricing/BlackScholes.hpp"
#include "../src/payoffs/VanillaPayoffs.hpp"

int main() {

    double S0 = 100;
    double K  = 100;
    double r  = 0.05;
    double sigma = 0.2;
    double T = 1.0;

    GeometricBrownianMotion gbm(r, sigma);

    auto paths =
        MonteCarloSimulator<EulerMaruyama>::simulate(
            gbm, S0, T, 0.01, 200000);

    CallPayoff payoff(K);

    double mc =
        EuropeanPricer::price(paths, payoff, r, T);

    double bs =
        BlackScholes::callPrice(S0, K, r, sigma, T);

    std::cout << "Monte Carlo price: " << mc << std::endl;
    std::cout << "Black-Scholes:    " << bs << std::endl;
}