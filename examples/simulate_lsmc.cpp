#include <iostream>
#include <iomanip>
#include <vector>
#include "../src/pricing/lsmc/AmericanPricer.hpp"
#include "../src/pricing/lsmc/BermudanPricer.hpp"
#include "../src/pricing/lsmc/EarlyExercisePremium.hpp"
#include "../src/payoffs/VanillaPayoffs.hpp"
#include "../src/pricing/BlackScholes.hpp"

static void separator() {
    std::cout << std::string(70, '-') << "\n";
}

int main() {
    std::cout << std::fixed << std::setprecision(4);

    const double spot    = 100.0;
    const double K       = 100.0;
    const double rate    = 0.05;
    const double T       = 1.0;
    const double flatVol = 0.20;
    const int    nSteps  = 50;    // exercise opportunities for American
    const int    nPaths  = 20000;

    // Rough Bergomi parameters (calibrated-style: H=0.12 as in literature)
    RoughBergomiModel::Parameters rbParams{
        /*xi0=*/  0.04,   // initial variance = 20% vol
        /*eta=*/  1.5,    // vol-of-vol
        /*rho=*/  -0.7,   // leverage
        /*hurst=*/0.12    // rough exponent
    };

    PutPayoff put(K);
    CallPayoff call(K);

    // ── LSMC config: augmented state for rough vol ────────────────────────
    AmericanPricer::Config roughCfg;
    roughCfg.degree          = 3;
    roughCfg.basis           = LSMCBasis::Laguerre;
    roughCfg.includeVol      = true;   // use v_t in regression
    roughCfg.includeVolterra = true;   // use Y_t in regression
    roughCfg.regularisation  = 1e-5;
    roughCfg.minITM          = 100;

    // Standard config: Markovian state only (S, t) — for comparison
    AmericanPricer::Config standardCfg;
    standardCfg.degree          = 3;
    standardCfg.basis           = LSMCBasis::Laguerre;
    standardCfg.includeVol      = false;
    standardCfg.includeVolterra = false;
    standardCfg.regularisation  = 1e-5;
    standardCfg.minITM          = 100;

    // ── 1. American put ───────────────────────────────────────────────────
    separator();
    std::cout << "American Put  S=" << spot << " K=" << K
              << " r=" << rate << " T=" << T << "\n";
    std::cout << "rBergomi: xi0=" << rbParams.xi0
              << " eta=" << rbParams.eta
              << " rho=" << rbParams.rho
              << " H=" << rbParams.hurst << "\n\n";

    auto americanRough    = AmericanPricer::price(rbParams, put, spot, rate, T,
                                nSteps, nPaths, roughCfg, 42);
    auto americanStandard = AmericanPricer::price(rbParams, put, spot, rate, T,
                                nSteps, nPaths, standardCfg, 42);
    const double bsEuro   = BlackScholes::putPrice(spot, K, rate, flatVol, T);

    std::cout << std::left << std::setw(40) << "LSMC (augmented: S,v,Y)"
              << std::right << std::setw(10) << americanRough.price
              << " ± " << americanRough.stdErr << "\n";
    std::cout << std::left << std::setw(40) << "LSMC (standard: S only)"
              << std::right << std::setw(10) << americanStandard.price
              << " ± " << americanStandard.stdErr << "\n";
    std::cout << std::left << std::setw(40) << "European MC (same paths)"
              << std::right << std::setw(10) << americanRough.europeanPrice << "\n";
    std::cout << std::left << std::setw(40) << "European BS (flat vol 20%)"
              << std::right << std::setw(10) << bsEuro << "\n";
    std::cout << "\n";
    std::cout << "Early exercise premium (rough, vs Euro MC):  "
              << americanRough.earlyExercisePremium << "\n";
    std::cout << "Augmented vs standard state impact:         "
              << (americanRough.price - americanStandard.price) << "\n";
    std::cout << "\nRegression R^2 at selected steps:\n";
    const int nDisplay = std::min(5, static_cast<int>(
        americanRough.rSquaredByStep.size()));
    for (int s = nSteps / nDisplay; s <= nSteps - 1;
             s += nSteps / nDisplay) {
        std::cout << "  step=" << std::setw(3) << s
                  << "  t=" << std::setprecision(3) << s * T / nSteps
                  << "  R2=" << std::setprecision(4)
                  << americanRough.rSquaredByStep[s] << "\n";
    }

    // ── 2. Bermudan put: quarterly exercise ───────────────────────────────
    separator();
    std::cout << "\nBermudan Put (quarterly: 4 exercise dates)\n\n";

    const std::vector<double> quarterlyDates = {0.25, 0.50, 0.75, 1.00};
    auto bermudanQ = BermudanPricer::price(
        rbParams, put, spot, rate, T,
        nSteps, nPaths, quarterlyDates, roughCfg, 42);

    std::cout << std::left << std::setw(40) << "Bermudan (quarterly, rough)"
              << std::right << std::setw(10) << bermudanQ.price
              << " ± " << bermudanQ.stdErr << "\n";
    std::cout << std::left << std::setw(40) << "American (50 dates, rough)"
              << std::right << std::setw(10) << americanRough.price << "\n";
    std::cout << std::left << std::setw(40) << "European"
              << std::right << std::setw(10) << bermudanQ.europeanPrice << "\n";
    std::cout << "\nPer-exercise-date diagnostics:\n";
    std::cout << std::left
              << std::setw(8) << "  Date"
              << std::setw(10) << "ExBdry"
              << std::setw(10) << "ExFrac"
              << std::setw(10) << "R^2" << "\n";
    for (const auto& d : bermudanQ.byDate) {
        std::cout << std::setw(8)  << d.time
                  << std::setw(10) << d.exerciseBoundary
                  << std::setw(10) << d.exerciseFraction
                  << std::setw(10) << d.rSquared << "\n";
    }

    // ── 3. Monthly Bermudan ───────────────────────────────────────────────
    separator();
    std::cout << "\nBermudan Put (monthly: 12 exercise dates)\n\n";

    std::vector<double> monthlyDates;
    for (int m = 1; m <= 12; ++m)
        monthlyDates.push_back(m / 12.0);

    auto bermudanM = BermudanPricer::price(
        rbParams, put, spot, rate, T,
        nSteps, nPaths, monthlyDates, roughCfg, 42);

    std::cout << std::left << std::setw(40) << "Bermudan (monthly, rough)"
              << std::right << std::setw(10) << bermudanM.price
              << " ± " << bermudanM.stdErr << "\n";

    // ── 4. Early exercise premium decomposition ───────────────────────────
    separator();
    std::cout << "\nEarly Exercise Premium Decomposition\n\n";

    auto decomp = EarlyExercisePremium::analyzePut(
        rbParams, spot, K, rate, flatVol, T,
        nSteps, nPaths, roughCfg, 42);

    std::cout << std::left << std::setw(40) << "American (rough Bergomi)"
              << std::right << decomp.americanPrice << "\n";
    std::cout << std::left << std::setw(40) << "European MC"
              << std::right << decomp.europeanMC << "\n";
    std::cout << std::left << std::setw(40) << "European BS (flat 20%)"
              << std::right << decomp.europeanBS << "\n";
    std::cout << std::left << std::setw(40) << "Barone-Adesi-Whaley (flat 20%)"
              << std::right << decomp.americanBAW << "\n";
    std::cout << "\n";
    std::cout << std::left << std::setw(40) << "Premium vs Euro MC"
              << std::right << decomp.premiumVsEuroMC << "\n";
    std::cout << std::left << std::setw(40) << "Premium vs Euro BS"
              << std::right << decomp.premiumVsEuroBS << "\n";
    std::cout << std::left << std::setw(40) << "vs BAW (model impact)"
              << std::right << decomp.premiumVsBAW << "\n";
    std::cout << std::left << std::setw(40) << "Rough vol impact on EEP"
              << std::right << decomp.roughVolImpact << "\n";

    // ── 5. Basis function comparison ──────────────────────────────────────
    separator();
    std::cout << "\nBasis Function Comparison (American Put, rough state)\n\n";

    for (auto [name, basisType] : std::vector<std::pair<std::string, LSMCBasis::Type>>{
            {"Laguerre",  LSMCBasis::Laguerre},
            {"Monomial",  LSMCBasis::Monomial},
            {"Chebyshev", LSMCBasis::Chebyshev}}) {
        AmericanPricer::Config c = roughCfg;
        c.basis = basisType;
        auto res = AmericanPricer::price(rbParams, put, spot, rate, T,
                                         nSteps, nPaths, c, 42);
        std::cout << std::left  << std::setw(16) << name
                  << std::right << std::setw(10) << res.price
                  << " ± " << res.stdErr << "\n";
    }

    // ── 6. Convergence in nPaths ──────────────────────────────────────────
    separator();
    std::cout << "\nConvergence in Number of Paths (American Put, rough Laguerre)\n\n";
    std::cout << std::left << std::setw(10) << "Paths"
              << std::right << std::setw(10) << "Price"
              << std::setw(12) << "StdErr" << "\n";

    for (int np : {2000, 5000, 10000, 20000, 50000}) {
        auto res = AmericanPricer::price(rbParams, put, spot, rate, T,
                                          nSteps, np, roughCfg, 42);
        std::cout << std::left << std::setw(10) << np
                  << std::right << std::setw(10) << res.price
                  << std::setw(12) << res.stdErr << "\n";
    }

    separator();
    return 0;
}