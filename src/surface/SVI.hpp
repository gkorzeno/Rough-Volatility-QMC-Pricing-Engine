// src/surface/SVI.hpp
#pragma once
#include <cmath>
#include <vector>
#include <stdexcept>
#include <functional>
#include <algorithm>
#include <array>

struct SVIParams {
    double a;      // vertical translation
    double b;      // slope/curvature
    double rho;    // rotation (skew)
    double m;      // horizontal translation
    double sigma;  // ATM curvature
};

class SVI {
public:
    template<typename T>
    static T clampValue(T x, T lo, T hi) {
        return std::max(lo, std::min(x, hi));
    }

    // Total variance at log-moneyness k = log(K/F)
    static double totalVariance(const SVIParams& p, double k) {
        double disc = k - p.m;
        return p.a + p.b * (p.rho * disc
               + std::sqrt(disc*disc + p.sigma*p.sigma));
    }

    // Implied vol at strike K, forward F, expiry T
    static double impliedVol(
        const SVIParams& p,
        double K, double F, double T)
    {
        double k = std::log(K / F);
        double w = totalVariance(p, k);
        if (w <= 0.0)
            throw std::runtime_error("SVI: negative total variance");
        return std::sqrt(w / T);
    }

    // Check butterfly arbitrage condition (g(k) >= 0)
    static bool isArbFree(const SVIParams& p, double k) {
        double disc  = k - p.m;
        double denom = std::sqrt(disc*disc + p.sigma*p.sigma);
        double w     = totalVariance(p, k);
        double dw    = p.b * (p.rho + disc / denom);
        double d2w   = p.b * p.sigma*p.sigma / (denom*denom*denom);

        double g = (1.0 - 0.5 * k * dw / w) * (1.0 - 0.5 * k * dw / w)
                 - 0.25 * dw*dw * (0.25 + 1.0/w)
                 + 0.5 * d2w;
        return g >= 0.0;
    }

    // Fit SVI to market vols via Nelder-Mead
    static SVIParams fit(
        const std::vector<double>& strikes,
        const std::vector<double>& marketVols,
        double F, double T,
        SVIParams init = {0.04, 0.1, -0.5, 0.0, 0.1})
    {
        // Objective: sum of squared vol errors
        auto objective = [&](const SVIParams& p) {
            double err = 0.0;
            for (size_t i = 0; i < strikes.size(); i++) {
                double modelVol = 0.0;
                try {
                    modelVol = impliedVol(p, strikes[i], F, T);
                } catch (...) {
                    return 1e12;
                }
                double diff     = modelVol - marketVols[i];
                err += diff * diff;
            }
            return err;
        };

        return nelderMead(objective, init);
    }

private:
    // Simple Nelder-Mead for 5-parameter SVI fit
    static SVIParams nelderMead(
        std::function<double(const SVIParams&)> f,
        SVIParams init,
        int maxIter = 200,
        double tol  = 1e-8)
    {
        // Represent SVIParams as 5D vector for simplex operations
        using Vec = std::array<double, 5>;

        auto toVec = [](const SVIParams& p) -> Vec {
            return {p.a, p.b, p.rho, p.m, p.sigma};
        };
        auto toParams = [](const Vec& v) -> SVIParams {
            // Clamp to valid SVI domain
            return {v[0],
                    std::max(v[1], 1e-6),
                    clampValue(v[2], -0.999, 0.999),
                    v[3],
                    std::max(v[4], 1e-6)};
        };
        auto fv = [&](const Vec& v) { return f(toParams(v)); };

        // Initialize simplex
        int n = 5;
        std::vector<Vec> simplex(n + 1);
        simplex[0] = toVec(init);
        for (int i = 0; i < n; i++) {
            simplex[i+1]    = simplex[0];
            simplex[i+1][i] *= (std::abs(simplex[0][i]) > 1e-6)
                               ? 1.05 : 0.00025;
        }

        for (int iter = 0; iter < maxIter; iter++) {
            // Sort by function value
            std::sort(simplex.begin(), simplex.end(),
                [&](const Vec& a, const Vec& b){ return fv(a) < fv(b); });

            if (fv(simplex.back()) - fv(simplex[0]) < tol)
                break;

            // Centroid (exclude worst)
            Vec cen = {};
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    cen[j] += simplex[i][j] / n;

            // Reflect
            Vec refl;
            for (int j = 0; j < n; j++)
                refl[j] = cen[j] + (cen[j] - simplex.back()[j]);

            if (fv(refl) < fv(simplex[0])) {
                // Expand
                Vec exp;
                for (int j = 0; j < n; j++)
                    exp[j] = cen[j] + 2.0*(cen[j] - simplex.back()[j]);
                simplex.back() = (fv(exp) < fv(refl)) ? exp : refl;
            } else if (fv(refl) < fv(simplex[n-1])) {
                simplex.back() = refl;
            } else {
                // Contract
                Vec con;
                for (int j = 0; j < n; j++)
                    con[j] = cen[j] + 0.5*(simplex.back()[j] - cen[j]);
                if (fv(con) < fv(simplex.back())) {
                    simplex.back() = con;
                } else {
                    // Shrink
                    for (int i = 1; i <= n; i++)
                        for (int j = 0; j < n; j++)
                            simplex[i][j] = simplex[0][j]
                                + 0.5*(simplex[i][j] - simplex[0][j]);
                }
            }
        }

        return toParams(simplex[0]);
    }
};
