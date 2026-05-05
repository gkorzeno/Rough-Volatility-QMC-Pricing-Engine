#include <iostream>
#include "../src/stochasticProcess/GeometricBrownianMotion.hpp"
#include "../src/simulators/PathSimulator.hpp"
#include "../src/integrators/EulerMaruyama.hpp"
#include "../src/strategies/ShortStraddleStrategy.hpp"
#include "../src/strategies/Backtester.hpp"

int main() {
    const double S0    = 100.0;
    const double K     = 100.0;
    const double r     = 0.05;
    const double sigma = 0.20;
    const double T     = 1.0;
    const double dt    = 0.01;
    const int    paths = 50000;

    GeometricBrownianMotion gbm(r, sigma);
    auto allPaths = PathSimulator<EulerMaruyama>::simulate(gbm, S0, T, dt, paths);

    // Price straddle at the true vol — sell it at implied == realised
    ShortStraddleStrategy strategy(K, sigma, /*tc=*/0.001);
    Backtester bt(r, T);

    const auto result = bt.run(strategy, allPaths);
    bt.report(strategy, result);

    // Now test vol mis-pricing: sell at 25% IV when realised is 20%
    ShortStraddleStrategy richStraddle(K, /*impliedVol=*/0.25, 0.001);
    const auto richResult = bt.run(richStraddle, allPaths);
    std::cout << "\n-- Sold at rich IV (25%) --\n";
    bt.report(richStraddle, richResult);

    return 0;
}