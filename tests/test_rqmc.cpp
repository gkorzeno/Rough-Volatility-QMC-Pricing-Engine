#include <iostream>
#include <iomanip>
#include "../src/stochasticProcess/GeometricBrownianMotion.hpp"
#include "../src/integrators/EulerMaruyama.hpp"
#include "../src/payoffs/VanillaPayoffs.hpp"
#include "../src/pricing/BlackScholes.hpp"
#include "../src/simulators/ParallelRQMCSimulator.hpp"

int main() {
    // Model params
    double S0 = 100.0;
    double K  = 100.0;
    double r  = 0.05;
    double sigma = 0.2;
    double T  = 1.0;

    double dt = 0.01;

    // RQMC params
    int pointsPerRep = 4096;
    int replicates   = 16;

    // Direction numbers file (Joe-Kuo table shipped in this repo)
    std::string dirFile = "docs/new-joe-kuo-6.21201.txt";

    // Process + payoff
    GeometricBrownianMotion gbm(r, sigma);
    CallPayoff payoff(K);

    std::cout << "=== Parallel RQMC Test ===\n\n";

    auto result = ParallelRQMCSimulator<EulerMaruyama>::priceEuropean(
        gbm,
        payoff,
        S0,
        r,
        T,
        dt,
        pointsPerRep,
        replicates,
        dirFile,
        true   // Brownian bridge ON
    );

    double bs = BlackScholes::callPrice(S0, K, r, sigma, T);

    std::cout << std::fixed << std::setprecision(6);

    std::cout << "RQMC Price      : " << result.price << "\n";
    std::cout << "Std Error       : " << result.stdErr << "\n";
    std::cout << "Black-Scholes   : " << bs << "\n";
    std::cout << "Absolute Error  : " << std::abs(result.price - bs) << "\n";

    std::cout << "\nDetails:\n";
    std::cout << "Replicates      : " << result.replicates << "\n";
    std::cout << "Points/Replicate: " << result.pointsPerReplicate << "\n";

    return 0;
}