#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>

#include "../src/greeks/AdjointGreeks.hpp"
#include "../src/pricing/BlackScholes.hpp"

int main() {
    std::cout << std::unitbuf;
    const double S0 = 100.0;
    const double K = 100.0;
    const double r = 0.05;
    const double sigma = 0.20;
    const double T = 1.0;
    const double dt = 0.01;
    const int paths = 100000;

    const Greeks aad = AdjointGreeks::computeCall(K, S0, r, sigma, T, dt, paths);
    const Greeks bs = BlackScholes::greeks(S0, K, r, sigma, T);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "=== Adjoint Greeks Test ===\n\n";
    std::cout << "AAD Price  : " << aad.price  << " | BS: " << bs.price  << "\n";
    std::cout << "AAD Delta  : " << aad.delta  << " | BS: " << bs.delta  << "\n";
    std::cout << "AAD Vega   : " << aad.vega   << " | BS: " << bs.vega   << "\n";
    std::cout << "AAD Rho    : " << aad.rho    << " | BS: " << bs.rho    << "\n";
    std::cout << "\n";
    std::cout << "|Price err| : " << std::abs(aad.price - bs.price) << "\n";
    std::cout << "|Delta err| : " << std::abs(aad.delta - bs.delta) << "\n";
    std::cout << "|Vega err|  : " << std::abs(aad.vega - bs.vega) << "\n";
    std::cout << "|Rho err|   : " << std::abs(aad.rho - bs.rho) << "\n";

    assert(std::abs(aad.price - bs.price) < 0.12);
    assert(std::abs(aad.delta - bs.delta) < 0.025);
    assert(std::abs(aad.vega - bs.vega) < 1.20);
    assert(std::abs(aad.rho - bs.rho) < 0.80);

    std::cout << "\nAdjoint first-order Greeks are within tolerance.\n";
    return 0;
}
