// examples/simulate_barrier.cpp
#include <iostream>
#include "../src/stochasticProcess/GeometricBrownianMotion.hpp"
#include "../src/barriers/Barrier.hpp"
#include "../src/payoffs/BarrierPayoff.hpp"
#include "../src/pricing/BarrierPricer.hpp"
#include "../src/pricing/BarrierAnalytical.hpp"
#include "../src/pricing/BlackScholes.hpp"
#include "../src/utils/Benchmark.hpp"

int main() {

    double S0    = 100.0;
    double K     = 100.0;
    double H_up  = 120.0;
    double H_dn  = 80.0;
    double r     = 0.05;
    double sigma = 0.2;
    double T     = 1.0;
    double dt    = 0.01;
    int    paths = 100000;
    // Check if random_device is working
    std::cout << "-- RNG diagnostics --\n";
    Random r1, r2, r3;
    std::cout << "Three normals from default RNG: "
          << r1.normal() << ", "
          << r2.normal() << ", "
          << r3.normal() << "\n";

    Random s1(42), s2(42);
    std::cout << "Two normals from seed=42 (should match): "
          << s1.normal() << ", " << s2.normal() << "\n";

    // Check magnitude — should rarely exceed 3
    Random test(123);
    int extreme = 0;
    for (int i = 0; i < 10000; i++)
        if (std::abs(test.normal()) > 3.0) extreme++;
    std::cout << "Draws > 3 sigma (expect ~27/10000): " << extreme << "\n\n";

    // Quick sanity check in simulate_barrier.cpp before the main runs
    std::cout << "-- Bridge probability sanity check --\n";
    // S=100, S=110, H=120, sigma=0.2, dt=0.01
    // Both below barrier — should be a small probability
    double p1 = BarrierCrossing::crossingProbability(100.0, 110.0, 120.0, 0.2, 0.01);
    // S=100, S=119, H=120 — very close, should be higher
    double p2 = BarrierCrossing::crossingProbability(100.0, 119.0, 120.0, 0.2, 0.01);
    // S=100, S=121 — definite cross, should be 1.0
    double p3 = BarrierCrossing::crossingProbability(100.0, 121.0, 120.0, 0.2, 0.01);
    std::cout << "p(100->110, H=120): " << p1 << " (expect ~0)\n";
    std::cout << "p(100->119, H=120): " << p2 << " (expect small but nonzero)\n";
    std::cout << "p(100->121, H=120): " << p3 << " (expect 1.0)\n";

    GeometricBrownianMotion gbm(r, sigma);
    BarrierCallPayoff call(K);

    double vanillaPrice = BlackScholes::callPrice(S0, K, r, sigma, T);
    std::cout << "Vanilla call (BS): " << vanillaPrice << "\n\n";

    // Replace the knock-out sanity check with this version
    std::cout << "-- Knock-out sanity check --\n";
    std::cout << "r=" << r << " sigma=" << sigma << " dt=" << dt
          << " T=" << T << " H=" << H_up << "\n";

    // Print first 5 steps of first path to verify GBM step size
    {
        Random rng(42);
        double S = S0;
        std::cout << "First 10 steps of path (expect small moves ~1% each):\n";
        for (int j = 0; j < 10; j++) {
            double z    = rng.normal();
            double Snew = S * std::exp((r - 0.5*sigma*sigma)*dt + sigma*std::sqrt(dt)*z);
            std::cout << "  j=" << j << " z=" << z 
                  << " S=" << S << " -> " << Snew 
                  << " (move=" << (Snew/S - 1)*100 << "%)\n";
            S = Snew;
        }
    }

    int crossCount = 0;
    {
        Random rng(42);
        int testPaths = 100000;
        int testSteps = static_cast<int>(T / dt);
        for (int i = 0; i < testPaths; i++) {
            double S = S0;
            bool hit = false;
            for (int j = 0; j < testSteps; j++) {
                double z = rng.normal();
                S = S * std::exp((r - 0.5*sigma*sigma)*dt + sigma*std::sqrt(dt)*z);
                if (S >= H_up) { hit = true; break; }
            }
            if (hit) crossCount++;
        }
        std::cout << "Discrete crossing rate: "
              << 100.0 * crossCount / testPaths << "%\n\n";
    }

    crossCount = 0;

    {
        Random rng(42);
        int testPaths = 100000;
        int testSteps = static_cast<int>(T / dt);
        for (int i = 0; i < testPaths; i++) {
            double S = S0;
            bool hit = false;
            for (int j = 0; j < testSteps; j++) {
                double z = rng.normal();
                S = S * std::exp((r - 0.5*sigma*sigma)*dt + sigma*std::sqrt(dt)*z);
                if (S >= H_up) { hit = true; break; }
            }
            if (hit) crossCount++;
        }
        std::cout << "Discrete crossing rate (no bridge): "
                << 100.0 * crossCount / testPaths << "%\n";
        std::cout << "Expected: ~17%\n\n";
    }

    // Test all four barrier types
    struct TestCase {
        std::string    name;
        Barrier        barrier;
        BarrierType    type;
    };

    std::vector<TestCase> cases = {
        {"Up-and-Out   (H=120)", Barrier(H_up, BarrierType::UpAndOut),   BarrierType::UpAndOut  },
        {"Up-and-In    (H=120)", Barrier(H_up, BarrierType::UpAndIn),    BarrierType::UpAndIn   },
        {"Down-and-Out (H=80)", Barrier(H_dn, BarrierType::DownAndOut), BarrierType::DownAndOut},
        {"Down-and-In  (H=80)", Barrier(H_dn, BarrierType::DownAndIn),  BarrierType::DownAndIn },
    };

    for (auto& tc : cases) {
        BarrierPricer::Result res_bridge, res_discrete;

        Benchmark::run(tc.name + " bridge  ", [&]() {
            res_bridge = BarrierPricer::price(
                gbm, call, tc.barrier, S0, r, sigma, T, dt, paths, true);
        });

        Benchmark::run(tc.name + " discrete", [&]() {
            res_discrete = BarrierPricer::price(
                gbm, call, tc.barrier, S0, r, sigma, T, dt, paths, false);
        });

        double analytical = BarrierAnalytical::price(
            S0, K, tc.barrier.H, r, sigma, T, tc.type);

        std::cout << "\n" << tc.name << "\n";
        std::cout << "  Analytical:          " << analytical          << "\n";
        std::cout << "  MC (bridge):         " << res_bridge.price
                  << " ± " << res_bridge.stdErr                       << "\n";
        std::cout << "  MC (discrete):       " << res_discrete.price
                  << " ± " << res_discrete.stdErr                     << "\n";
        std::cout << "  Knocked out (bridge):" << res_bridge.fractionKnockedOut * 100
                  << "%\n";
    }

    // Key identity: Up-and-In + Up-and-Out = Vanilla (same H)
    std::cout << "\n-- Parity check (Up-and-In + Up-and-Out = Vanilla) --\n";
    auto uao = BarrierPricer::price(gbm, call,
        Barrier(H_up, BarrierType::UpAndOut), S0, r, sigma, T, dt, paths, true);
    auto uai = BarrierPricer::price(gbm, call,
        Barrier(H_up, BarrierType::UpAndIn),  S0, r, sigma, T, dt, paths, true);
    std::cout << "  Up-and-Out + Up-and-In = " << uao.price + uai.price << "\n";
    std::cout << "  Vanilla                = " << vanillaPrice           << "\n";

    return 0;
}