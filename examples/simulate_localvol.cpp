// examples/simulate_localvol.cpp
#include <iostream>
#include <iomanip>
#include <vector>
#include "../src/surface/VolSurface.hpp"
#include "../src/surface/LocalVolSurface.hpp"
#include "../src/surface/SVI.hpp"
#include "../src/stochasticProcess/CEV.hpp"
#include "../src/stochasticProcess/LocalVolProcess.hpp"
#include "../src/calibration/ImpliedVolSolver.hpp"
#include "../src/simulators/CEVSimulator.hpp"
#include "../src/simulators/ParallelMonteCarloSimulator.hpp"
#include "../src/payoffs/VanillaPayoffs.hpp"
#include "../src/pricing/EuropeanPricer.hpp"
#include "../src/pricing/BlackScholes.hpp"
#include "../src/integrators/EulerMaruyama.hpp"

int main() {
    double S0    = 100.0;
    double r     = 0.05;
    double T     = 1.0;
    double dt    = 0.01;
    int    paths = 50000;

    // ── 1. CEV model ─────────────────────────────────────────────────────
    std::cout << "=== CEV Model ===\n";
    std::cout << std::fixed << std::setprecision(4);
    std::cout << std::setw(8)  << "beta"
              << std::setw(12) << "ATM price"
              << std::setw(12) << "OTM price"
              << std::setw(12) << "BS (beta=1)\n";
    std::cout << std::string(44, '-') << "\n";

    double bsATM = BlackScholes::callPrice(S0, 100.0, r, 0.2, T);
    double bsOTM = BlackScholes::callPrice(S0, 110.0, r, 0.2, T);

    for (double beta : {0.3, 0.5, 0.7, 1.0}) {
        CEV cev(r, 0.2, beta);
        CallPayoff atm(100.0), otm(110.0);

        auto res = CEVSimulator::simulate(cev, S0, T, dt, paths);

        auto [atmPrice, atmErr] = EuropeanPricer::priceWithError(res, atm, r, T);
        auto [otmPrice, otmErr] = EuropeanPricer::priceWithError(res, otm, r, T);

        std::cout << std::setw(8)  << beta
                  << std::setw(12) << atmPrice
                  << std::setw(12) << otmPrice
                  << std::setw(12) << (beta == 1.0 ? bsATM : 0.0) << "\n";
    }

    // ── 2. SVI smile fitting ─────────────────────────────────────────────
    std::cout << "\n=== SVI Smile Fitting ===\n";

    double F = S0 * std::exp(r * T);
    std::vector<double> strikes     = {80, 85, 90, 95, 100, 105, 110, 115, 120};
    std::vector<double> marketVols;

    for (double K : strikes) {
        double k   = std::log(K / F);
        double vol = 0.2 - 0.1*k + 0.05*k*k;
        marketVols.push_back(std::max(vol, 0.01));
    }

    SVIParams fitted = SVI::fit(strikes, marketVols, F, T);

    std::cout << "Fitted SVI params:\n";
    std::cout << "  a="     << fitted.a
              << " b="      << fitted.b
              << " rho="    << fitted.rho
              << " m="      << fitted.m
              << " sigma="  << fitted.sigma << "\n\n";

    std::cout << std::setw(8)  << "Strike"
              << std::setw(12) << "Market Vol"
              << std::setw(12) << "SVI Vol"
              << std::setw(10) << "Error\n";
    std::cout << std::string(42, '-') << "\n";

    for (size_t i = 0; i < strikes.size(); i++) {
        double sviVol = SVI::impliedVol(fitted, strikes[i], F, T);
        std::cout << std::setw(8)  << strikes[i]
                  << std::setw(12) << marketVols[i]
                  << std::setw(12) << sviVol
                  << std::setw(10) << (sviVol - marketVols[i]) << "\n";
    }

    // Arbitrage check across strikes
    std::cout << "\nButterfly arbitrage check (g >= 0 required):\n";
    bool arbFree = true;
    for (double K : strikes) {
        double k = std::log(K / F);
        bool ok  = SVI::isArbFree(fitted, k);
        if (!ok) arbFree = false;
        std::cout << "  K=" << K << " k=" << k
                  << " arb-free=" << (ok ? "yes" : "NO") << "\n";
    }
    std::cout << "Overall: " << (arbFree ? "arbitrage-free\n" : "ARBITRAGE DETECTED\n");

    // ── 3. Vol surface + Dupire local vol ────────────────────────────────
    std::cout << "\n=== Local Vol Surface ===\n";

    std::vector<double> expiries     = {0.25, 0.5, 1.0, 2.0};
    std::vector<std::vector<double>> surfaceVols;

    // IMPORTANT: SVI is parameterized in TOTAL variance w(k,T), not sigma.
    // Re-using one raw SVI parameter set across all maturities makes
    // implied vol scale as sqrt(w/T), which explodes at short T.
    // Build a realistic implied-vol surface directly by strike/time instead.
    for (double t : expiries) {
        double Ft = S0 * std::exp(r * t);
        std::vector<double> row;
        for (double K : strikes) {
            double k = std::log(K / Ft);
            // Base smile around 20% with mild skew + weak term-structure.
            double baseVol = 0.20 - 0.10 * k + 0.05 * k * k;
            double termAdj = 0.015 * std::exp(-t); // slightly higher short-end
            row.push_back(std::max(0.05, baseVol + termAdj));
        }
        surfaceVols.push_back(row);
    }

    VolSurface      volSurf(strikes, expiries, surfaceVols);
    LocalVolSurface lvSurf(volSurf, S0, r);

    // Diagnostics — check surface and Dupire values before simulating
    std::cout << "\n-- Implied vol surface check --\n";
    std::cout << std::setw(8)  << "K"
              << std::setw(8)  << "T"
              << std::setw(14) << "Implied Vol"
              << std::setw(14) << "Local Vol\n";
    std::cout << std::string(44, '-') << "\n";
    for (double K : {90.0, 100.0, 110.0}) {
        for (double t : {0.25, 0.5, 1.0}) {
            double iv = volSurf.impliedVol(t, K);
            double lv = lvSurf.localVol(K, t);
            std::cout << std::setw(8)  << K
                      << std::setw(8)  << t
                      << std::setw(14) << iv
                      << std::setw(14) << lv << "\n";
        }
    }

    // Simulate local vol paths and price options
    std::cout << "\n-- Local vol option prices --\n";
    LocalVolProcess lvProcess(lvSurf, r);

    auto lvRes = ParallelMonteCarloSimulator<EulerMaruyama>::simulate(
        lvProcess, S0, T, dt, paths);

    std::cout << std::setw(10) << "Strike"
              << std::setw(14) << "Local Vol MC"
              << std::setw(14) << "BS (flat)"
              << std::setw(14) << "Difference\n";
    std::cout << std::string(52, '-') << "\n";

    for (double K : {90.0, 95.0, 100.0, 105.0, 110.0}) {
        CallPayoff call(K);
        auto [lvPrice, lvErr] = EuropeanPricer::priceWithError(lvRes, call, r, T);
        double bsPrice        = BlackScholes::callPrice(S0, K, r, 0.2, T);
        std::cout << std::setw(10) << K
                  << std::setw(14) << lvPrice
                  << std::setw(14) << bsPrice
                  << std::setw(14) << (lvPrice - bsPrice) << "\n";
    }

    // ── 4. Implied vol inversion ─────────────────────────────────────────
    std::cout << "\n=== Implied Vol Inversion ===\n";
    std::cout << std::setw(8)  << "Strike"
              << std::setw(12) << "BS Price"
              << std::setw(12) << "Impl Vol"
              << std::setw(12) << "True Vol"
              << std::setw(10) << "Error\n";
    std::cout << std::string(54, '-') << "\n";

    for (double K : {90.0, 95.0, 100.0, 105.0, 110.0}) {
        double trueVol = 0.2;
        double price   = BlackScholes::callPrice(S0, K, r, trueVol, T);
        double implVol = ImpliedVolSolver::solve(price, S0, K, r, T);
        std::cout << std::setw(8)  << K
                  << std::setw(12) << price
                  << std::setw(12) << implVol
                  << std::setw(12) << trueVol
                  << std::setw(10) << (implVol - trueVol) << "\n";
    }

    // ── 5. Back out implied vols from local vol MC prices ────────────────
    std::cout << "\n=== Implied Vols from Local Vol MC ===\n";
    std::cout << std::setw(8)  << "Strike"
              << std::setw(14) << "LV MC Price"
              << std::setw(14) << "Impl Vol"
              << std::setw(14) << "Flat Vol\n";
    std::cout << std::string(50, '-') << "\n";

    for (double K : {90.0, 95.0, 100.0, 105.0, 110.0}) {
        CallPayoff call(K);
        auto [lvPrice, lvErr] = EuropeanPricer::priceWithError(lvRes, call, r, T);

        double implVol = 0.0;
        try {
            implVol = ImpliedVolSolver::solve(lvPrice, S0, K, r, T);
        } catch (...) {
            implVol = -1.0;  // mark as failed
        }

        std::cout << std::setw(8)  << K
                  << std::setw(14) << lvPrice
                  << std::setw(14) << implVol
                  << std::setw(14) << 0.2 << "\n";
    }

    return 0;
}