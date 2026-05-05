#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include "../../surface/VolSurface.hpp"
#include "../../surface/SVI.hpp"

class SkewAnalyzer {
public:
    struct SkewMetrics {
        double expiry;
        double atmVol;
        double skew;          // (25d put IV - 25d call IV) / ATM
        double convexity;     // (25d put IV + 25d call IV - 2*ATM) / ATM
        double riskReversalIV; // put25 - call25
        double butterflyIV;   // (put25 + call25)/2 - atm
        SVIParams sviParams;
    };

    struct TermStructureMetrics {
        double shortExpiry, longExpiry;
        double shortATM,    longATM;
        double termSpread;    // longATM - shortATM
        double forwardVol;    // implied forward vol between the two expiries
        bool   termInversion; // short vol > long vol (inverted term structure)
    };

    struct Signal {
        enum Type { None, CalendarSpread, RiskReversal, Butterfly } type = None;
        double longExpiry, longStrike, longIV;
        double shortExpiry, shortStrike, shortIV;
        double expectedEdge;
        std::string description;
    };

    static SkewMetrics computeSkew(
        const VolSurface& surface,
        double expiry,
        double spot,
        double rate = 0.05)
    {
        const double F      = spot * std::exp(rate * expiry);
        const double atmVol = surface.impliedVol(expiry, spot);
        const double sqrtT  = std::sqrt(expiry);

        // Approximate 25-delta strikes via BS inverse
        // For a call with delta=0.25: K_c25 = F * exp(-d1 * sigma * sqrtT + ...)
        // Quick approximation: K_25c ≈ F * exp(+0.674 * sigma * sqrtT)
        //                      K_25p ≈ F * exp(-0.674 * sigma * sqrtT)
        const double bump = 0.674 * atmVol * sqrtT;
        const double K25c = F * std::exp( bump);
        const double K25p = F * std::exp(-bump);

        const double iv25c = surface.impliedVol(expiry, K25c);
        const double iv25p = surface.impliedVol(expiry, K25p);

        // Fit SVI to get smooth parameters
        const auto& strikes = surface.getStrikes();
        std::vector<double> vols;
        vols.reserve(strikes.size());
        for (double K : strikes)
            vols.push_back(surface.impliedVol(expiry, K));

        SVIParams svi = SVI::fit(strikes, vols, F, expiry);

        return {expiry, atmVol,
                (iv25p - iv25c) / atmVol,
                (iv25p + iv25c - 2.0 * atmVol) / atmVol,
                iv25p - iv25c,
                (iv25p + iv25c) / 2.0 - atmVol,
                svi};
    }

    static TermStructureMetrics computeTermStructure(
        const VolSurface& surface,
        double shortExpiry,
        double longExpiry,
        double spot)
    {
        const double shortATM = surface.impliedVol(shortExpiry, spot);
        const double longATM  = surface.impliedVol(longExpiry,  spot);

        // Forward vol: sigma_fwd^2 * (T2 - T1) = sigma_long^2 * T2 - sigma_short^2 * T1
        const double fwdVar = std::max(0.0,
            longATM  * longATM  * longExpiry
          - shortATM * shortATM * shortExpiry);
        const double fwdVol = (longExpiry > shortExpiry)
            ? std::sqrt(fwdVar / (longExpiry - shortExpiry))
            : longATM;

        return {shortExpiry, longExpiry,
                shortATM, longATM,
                longATM - shortATM,
                fwdVol,
                shortATM > longATM};
    }

    // Generate trading signals from surface analysis
    static std::vector<Signal> generateSignals(
        const VolSurface& surface,
        double spot, double rate,
        double skewThreshold    = 0.10,  // 10% normalised skew to trade RR
        double termSpreadThresh = 0.02)  // 2 vol points to trade calendar
    {
        const auto& expiries = surface.getExpiries();
        const auto& strikes  = surface.getStrikes();
        std::vector<Signal> signals;

        // ── Calendar spread signals ───────────────────────────────────────
        for (std::size_t i = 0; i + 1 < expiries.size(); ++i) {
            const auto ts = computeTermStructure(
                surface, expiries[i], expiries[i+1], spot);

            if (ts.termInversion || ts.termSpread > termSpreadThresh) {
                Signal sig;
                sig.type        = Signal::CalendarSpread;
                sig.longExpiry  = ts.termInversion ? ts.longExpiry  : ts.shortExpiry;
                sig.shortExpiry = ts.termInversion ? ts.shortExpiry : ts.longExpiry;
                sig.longStrike  = spot;
                sig.shortStrike = spot;
                sig.longIV      = ts.termInversion ? ts.longATM  : ts.shortATM;
                sig.shortIV     = ts.termInversion ? ts.shortATM : ts.longATM;
                sig.expectedEdge = std::abs(ts.termSpread);
                sig.description  = ts.termInversion
                    ? "Sell front/buy back (inverted term structure)"
                    : "Buy front/sell back (steep term structure)";
                signals.push_back(sig);
            }
        }

        // ── Risk reversal signals (skew too steep/flat) ───────────────────
        for (double T : expiries) {
            const auto sk = computeSkew(surface, T, spot, rate);
            if (std::abs(sk.skew) > skewThreshold) {
                const double F    = spot * std::exp(rate * T);
                const double bump = 0.674 * sk.atmVol * std::sqrt(T);

                Signal sig;
                sig.type        = Signal::RiskReversal;
                sig.longExpiry  = sig.shortExpiry = T;
                sig.expectedEdge = std::abs(sk.riskReversalIV);

                if (sk.skew > 0) {
                    // Put skew expensive: sell put25, buy call25
                    sig.shortStrike = F * std::exp(-bump);
                    sig.longStrike  = F * std::exp(+bump);
                    sig.shortIV     = surface.impliedVol(T, sig.shortStrike);
                    sig.longIV      = surface.impliedVol(T, sig.longStrike);
                    sig.description = "Sell put skew / buy call (skew too steep)";
                } else {
                    // Call skew expensive: sell call25, buy put25
                    sig.shortStrike = F * std::exp(+bump);
                    sig.longStrike  = F * std::exp(-bump);
                    sig.shortIV     = surface.impliedVol(T, sig.shortStrike);
                    sig.longIV      = surface.impliedVol(T, sig.longStrike);
                    sig.description = "Sell call skew / buy put (flat skew)";
                }
                signals.push_back(sig);
            }
        }

        // Sort by expected edge descending
        std::sort(signals.begin(), signals.end(),
            [](const Signal& a, const Signal& b) {
                return a.expectedEdge > b.expectedEdge; });

        return signals;
    }
};