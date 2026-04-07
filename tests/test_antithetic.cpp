#include <iostream>
#include <numeric>
#include <cmath>
#include <vector>

#include "../src/stochasticProcess/GeometricBrownianMotion.hpp"
#include "../src/simulators/MonteCarloSimulator.hpp"
#include "../src/simulators/AntitheticSimulator.hpp"
#include "../src/integrators/EulerMaruyama.hpp"
#include "../src/payoffs/VanillaPayoffs.hpp"

double variance(const std::vector<double>& x)
{
    double mean = std::accumulate(x.begin(), x.end(), 0.0) / x.size();

    double var = 0;
    for (double v : x)
        var += (v - mean)*(v - mean);

    return var / x.size();
}

double mean(const std::vector<double>& x)
{
    return std::accumulate(x.begin(), x.end(), 0.0) / x.size();
}

int main()
{
    GeometricBrownianMotion gbm(0.05, 0.2);

    double S0 = 100;
    double K  = 100;
    double T  = 1.0;
    double dt = 0.01;
    int paths = 100000;

    CallPayoff payoff(K);

    // simulate terminal prices
    auto paths_regular =
        MonteCarloSimulator<EulerMaruyama>::simulate(
            gbm, S0, T, dt, paths);

    auto paths_anti =
        AntitheticSimulator<EulerMaruyama>::simulate(
            gbm, S0, T, dt, paths);

    std::vector<double> reg_payoffs;
    std::vector<double> anti_payoffs;

    for(double s : paths_regular)
        reg_payoffs.push_back(payoff(s));

    // for(double s : paths_anti)
    //     anti_payoffs.push_back(payoff(s));

    for(int i = 0; i < paths_anti.size(); i += 2)
    {
        double p1 = payoff(paths_anti[i]);
        double p2 = payoff(paths_anti[i+1]);

        anti_payoffs.push_back(0.5 * (p1 + p2));
    }

    std::cout << "Regular payoff variance:    "
              << variance(reg_payoffs) << std::endl;

    std::cout << "Antithetic payoff variance: "
              << variance(anti_payoffs) << std::endl;

    std::cout << "Regular mean: "
              << mean(reg_payoffs) << std::endl;

    std::cout << "Antithetic mean: "
              << mean(anti_payoffs) << std::endl;
}