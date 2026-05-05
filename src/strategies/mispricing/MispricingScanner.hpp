#pragma once
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include "../../surface/VolSurface.hpp"
#include "../../surface/SVI.hpp"
#include "../../roughvol/RoughBergomiModel.hpp"
#include "../../surface/LocalVolSurface.hpp"
#include "../../calibration/ImpliedVolSolver.hpp"
#include "../../pricing/BlackScholes.hpp"

class MispricingScanner {
public:
    // How much each model disagrees with the market at a given (T, K) point
    struct MispricingPoint {
        double expiry;
        double strike;
        double marketIV;
        double sviIV;
        double rbergomiIV;
        double localvolIV;
        double sviMispricing;       // sviIV - marketIV
        double rbergomiMispricing;  // rbergomiIV - marketIV
        double localvolMispricing;  // localvolIV - marketIV
        double consensusMispricing; // mean of all model mispricings
    };

    struct RelativeValueTrade {
        double longExpiry,  longStrike;   // buy the undervalued leg
        double shortExpiry, shortStrike;  // sell the overvalued leg
        double longIV,      shortIV;
        double netEdge;     // shortIV - longIV (positive = positive carry)
        std::string description;
    };

    // Scan a vol surface for points where models consistently disagree with market
    static std::vector<MispricingPoint> scan(
        const VolSurface& market,
        const VolSurface& sviSurface,
        const VolSurface& rbergomiSurface,
        const VolSurface& localvolSurface,
        double minAbsMispricing = 0.005)   // 0.5 vol points threshold
    {
        const auto& expiries = market.getExpiries();
        const auto& strikes  = market.getStrikes();

        std::vector<MispricingPoint> results;

        for (double T : expiries) {
            for (double K : strikes) {
                const double mktIV  = market.impliedVol(T, K);
                const double sviIV  = sviSurface.impliedVol(T, K);
                const double rbIV   = rbergomiSurface.impliedVol(T, K);
                const double lvIV   = localvolSurface.impliedVol(T, K);

                const double sviMisp = sviIV  - mktIV;
                const double rbMisp  = rbIV   - mktIV;
                const double lvMisp  = lvIV   - mktIV;
                const double consensus = (sviMisp + rbMisp + lvMisp) / 3.0;

                if (std::abs(consensus) >= minAbsMispricing) {
                    results.push_back({T, K, mktIV, sviIV, rbIV, lvIV,
                                       sviMisp, rbMisp, lvMisp, consensus});
                }
            }
        }

        // Sort by absolute consensus mispricing descending
        std::sort(results.begin(), results.end(),
            [](const MispricingPoint& a, const MispricingPoint& b) {
                return std::abs(a.consensusMispricing) > std::abs(b.consensusMispricing);
            });

        return results;
    }

    // Replace the buildRVTrades method:
    static std::vector<RelativeValueTrade> buildRVTrades(
        const std::vector<MispricingPoint>& mispricings,
        int maxTrades = 5)
    {
        // cheap = market IV is BELOW model consensus (models say it's worth more)
        //         → buy these options (long the undervalued leg)
        // rich  = market IV is ABOVE model consensus (models say it's worth less)
        //         → sell these options (short the overvalued leg)
        std::vector<const MispricingPoint*> cheap, rich;
        for (const auto& mp : mispricings) {
        // consensusMispricing = (svi+rb+lv)/3 - market
        // positive  → models above market → market is cheap → BUY
        // negative  → models below market → market is rich  → SELL
            if (mp.consensusMispricing > 0)
                cheap.push_back(&mp);
            else
                rich.push_back(&mp);
        }

        // Sort cheap ascending by mispricing (least undervalued first — most
        // conservative long legs) and rich descending (most overvalued first)
        std::sort(cheap.begin(), cheap.end(),
            [](const MispricingPoint* a, const MispricingPoint* b) {
                return a->consensusMispricing < b->consensusMispricing; });
        std::sort(rich.begin(), rich.end(),
            [](const MispricingPoint* a, const MispricingPoint* b) {
                return a->consensusMispricing > b->consensusMispricing; });

        std::vector<RelativeValueTrade> trades;
        const int n = std::min({maxTrades,
                            static_cast<int>(cheap.size()),
                            static_cast<int>(rich.size())});

        for (int i = 0; i < n; ++i) {
            const auto* longLeg  = cheap[i];
            const auto* shortLeg = rich[i];

            trades.push_back({
                longLeg->expiry,   longLeg->strike,
                shortLeg->expiry,  shortLeg->strike,
                longLeg->marketIV, shortLeg->marketIV,
                shortLeg->marketIV - longLeg->marketIV,
                "long T=" + fmt(longLeg->expiry)  + " K=" + fmt(longLeg->strike)
                + " / short T=" + fmt(shortLeg->expiry) + " K=" + fmt(shortLeg->strike)
            });
        }
        return trades;
    }

private:
    static std::string fmt(double x) {
        char buf[32];
        std::snprintf(buf, sizeof(buf), "%.2f", x);
        return buf;
    }
};