// src/surface/VolSurface.hpp
#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>

// Stores implied vol as a grid over (T, K)
// Supports bilinear interpolation
class VolSurface {
private:
    template<typename T>
    static T clampValue(T x, T lo, T hi) {
        return std::max(lo, std::min(x, hi));
    }

    std::vector<double> strikes;   // K grid
    std::vector<double> expiries;  // T grid
    // vols[i][j] = implied vol at expiry i, strike j
    std::vector<std::vector<double>> vols;

public:
    VolSurface(
        const std::vector<double>& strikes_,
        const std::vector<double>& expiries_,
        const std::vector<std::vector<double>>& vols_)
        : strikes(strikes_), expiries(expiries_), vols(vols_)
    {
        if (vols_.size() != expiries_.size())
            throw std::runtime_error("VolSurface: expiry/vol dimension mismatch");
        for (auto& row : vols_)
            if (row.size() != strikes_.size())
                throw std::runtime_error("VolSurface: strike/vol dimension mismatch");
    }

    // Bilinear interpolation
    double impliedVol(double T, double K) const {
        // Clamp to grid boundaries
        T = clampValue(T, expiries.front(), expiries.back());
        K = clampValue(K, strikes.front(),  strikes.back());

        // Find bracketing indices
        int ti = upperIdx(expiries, T);
        int ki = upperIdx(strikes,  K);

        // Bilinear weights
        double t0 = expiries[ti-1], t1 = expiries[ti];
        double k0 = strikes[ki-1],  k1 = strikes[ki];
        double wt  = (T - t0) / (t1 - t0);
        double wk  = (K - k0) / (k1 - k0);

        return (1-wt)*(1-wk)*vols[ti-1][ki-1]
             + (1-wt)*  wk  *vols[ti-1][ki  ]
             +    wt *(1-wk)*vols[ti  ][ki-1]
             +    wt *   wk *vols[ti  ][ki  ];
    }

    // Total variance w(T,K) = sigma^2 * T — needed for Dupire
    double totalVariance(double T, double K) const {
        double sig = impliedVol(T, K);
        return sig * sig * T;
    }

    const std::vector<double>& getStrikes()  const { return strikes;  }
    const std::vector<double>& getExpiries() const { return expiries; }

private:
    static int upperIdx(const std::vector<double>& v, double x) {
        auto it = std::upper_bound(v.begin(), v.end(), x);
        int  i  = std::distance(v.begin(), it);
        return clampValue(i, 1, (int)v.size() - 1);
    }
};
