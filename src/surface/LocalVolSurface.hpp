// src/surface/LocalVolSurface.hpp
#pragma once
#include <cmath>
#include <vector>
#include "SVI.hpp"
#include "VolSurface.hpp"
#include "../pricing/BlackScholes.hpp"

class LocalVolSurface {
private:
    template<typename T>
    static T clampValue(T x, T lo, T hi) {
        return std::max(lo, std::min(x, hi));
    }

    VolSurface rawSurface;
    VolSurface smoothedSurface;
    double S0, r;
    double dK, dT;

    static double safeCallPrice(double S, double K, double r, double sigma, double T) {
        const double tEff = std::max(T, 1e-8);
        const double sigEff = clampValue(sigma, 1e-4, 3.0);
        return BlackScholes::callPrice(S, K, r, sigEff, tEff);
    }

    static std::vector<std::vector<double>> buildSviSmoothedVols(
        const VolSurface& surface,
        double S0,
        double r)
    {
        const std::vector<double>& strikes = surface.getStrikes();
        const std::vector<double>& expiries = surface.getExpiries();
        std::vector<std::vector<double>> smooth(
            expiries.size(),
            std::vector<double>(strikes.size(), 0.0));

        for (std::size_t ti = 0; ti < expiries.size(); ++ti) {
            const double T = expiries[ti];
            const double F = S0 * std::exp(r * T);
            std::vector<double> row(strikes.size(), 0.0);
            for (std::size_t ki = 0; ki < strikes.size(); ++ki)
                row[ki] = surface.impliedVol(T, strikes[ki]);

            try {
                const SVIParams fit = SVI::fit(strikes, row, F, T);
                for (std::size_t ki = 0; ki < strikes.size(); ++ki) {
                    double sviVol = row[ki];
                    try {
                        sviVol = SVI::impliedVol(fit, strikes[ki], F, T);
                    } catch (...) {
                        sviVol = row[ki];
                    }
                    smooth[ti][ki] = clampValue(sviVol, 0.01, 3.0);
                }
            } catch (...) {
                smooth[ti] = row;
            }
        }

        return smooth;
    }

    double smoothedImpliedVol(double T, double K) const {
        return smoothedSurface.impliedVol(T, K);
    }

    double callPriceFromSurface(double K, double T) const {
        const double kEff = clampValue(K, rawSurface.getStrikes().front(), rawSurface.getStrikes().back());
        const double tEff = clampValue(std::max(T, 1e-8), rawSurface.getExpiries().front(), rawSurface.getExpiries().back());
        const double sig = smoothedImpliedVol(tEff, kEff);
        return safeCallPrice(S0, kEff, r, sig, tEff);
    }

    double strikeSecondDerivative(double K, double T, double h) const {
        const double kMin = rawSurface.getStrikes().front();
        const double kMax = rawSurface.getStrikes().back();
        const double left = std::max(kMin, K - h);
        const double right = std::min(kMax, K + h);

        if (right - left < 1e-10)
            return 0.0;

        if (std::abs((K - left) - (right - K)) < 1e-10) {
            const double c0 = callPriceFromSurface(K, T);
            const double cUp = callPriceFromSurface(right, T);
            const double cDown = callPriceFromSurface(left, T);
            const double hh = right - K;
            return (cUp - 2.0 * c0 + cDown) / (hh * hh);
        }

        const double cL = callPriceFromSurface(left, T);
        const double c0 = callPriceFromSurface(K, T);
        const double cR = callPriceFromSurface(right, T);
        const double h1 = K - left;
        const double h2 = right - K;
        return 2.0 * (cL / (h1 * (h1 + h2))
                    - c0 / (h1 * h2)
                    + cR / (h2 * (h1 + h2)));
    }

public:
    LocalVolSurface(
        const VolSurface& surface_,
        double S0_, double r_,
        double dK_ = 0.0, double dT_ = 0.0)
        : rawSurface(surface_),
          smoothedSurface(
              surface_.getStrikes(),
              surface_.getExpiries(),
              buildSviSmoothedVols(surface_, S0_, r_)),
          S0(S0_), r(r_)
    {
        const std::vector<double>& K = rawSurface.getStrikes();
        const std::vector<double>& T = rawSurface.getExpiries();
        const double strikeRange = K.back() - K.front();
        const double expiryRange = T.back() - T.front();
        dK = (dK_ > 0) ? dK_ : std::max(0.5, 0.03 * strikeRange);
        dT = (dT_ > 0) ? dT_ : std::max(1.0 / 365.0, 0.08 * expiryRange);
    }

    double localVol(double K, double T) const {
        const double kMin = rawSurface.getStrikes().front();
        const double kMax = rawSurface.getStrikes().back();
        const double tMin = rawSurface.getExpiries().front();
        const double tMax = rawSurface.getExpiries().back();
        const double kEff = clampValue(K, kMin, kMax);
        const double tEff = clampValue(std::max(T, tMin), tMin, tMax);
        const double impliedFallback = smoothedImpliedVol(tEff, kEff);

        const double baseHk = std::max(0.5, std::min(dK, 0.15 * std::max(kEff, 1.0)));
        const double leftRoom = kEff - kMin;
        const double rightRoom = kMax - kEff;
        double hk = baseHk;
        if (leftRoom < hk || rightRoom < hk)
            hk = std::max(0.25, std::min(baseHk, std::max(leftRoom, rightRoom)));

        const double htForward = std::min(dT, std::max(1e-6, tMax - tEff));
        const double htBackward = std::min(dT, std::max(1e-6, tEff - tMin));

        const double c0 = callPriceFromSurface(kEff, tEff);
        double dCdT = 0.0;
        if (htForward > 1e-8 && htBackward > 1e-8) {
            const double cUpT = callPriceFromSurface(kEff, tEff + htForward);
            const double cDownT = callPriceFromSurface(kEff, tEff - htBackward);
            dCdT = (cUpT - cDownT) / (htForward + htBackward);
        } else if (htForward > 1e-8) {
            const double cUpT = callPriceFromSurface(kEff, tEff + htForward);
            dCdT = (cUpT - c0) / htForward;
        } else if (htBackward > 1e-8) {
            const double cDownT = callPriceFromSurface(kEff, tEff - htBackward);
            dCdT = (c0 - cDownT) / htBackward;
        }

        double dCdK = 0.0;
        const double left = std::max(kMin, kEff - hk);
        const double right = std::min(kMax, kEff + hk);
        if (right - left > 1e-10) {
            const double cUpK = callPriceFromSurface(right, tEff);
            const double cDownK = callPriceFromSurface(left, tEff);
            dCdK = (cUpK - cDownK) / (right - left);
        }
        const double d2CdK2 = strikeSecondDerivative(kEff, tEff, hk);

        const double numerator = dCdT + r * kEff * dCdK;
        const double denominator = 0.5 * kEff * kEff * d2CdK2;

        if (!std::isfinite(numerator) || !std::isfinite(denominator) || std::abs(denominator) < 1e-8)
            return impliedFallback;

        const double localVariance = numerator / denominator;
        if (!std::isfinite(localVariance) || localVariance <= 0.0)
            return impliedFallback;

        const double sigLoc = clampValue(std::sqrt(localVariance), 0.01, 2.0);
        return std::isfinite(sigLoc) ? sigLoc : impliedFallback;
    }

    double smoothedImpliedVolAt(double T, double K) const {
        return smoothedImpliedVol(T, K);
    }
};
