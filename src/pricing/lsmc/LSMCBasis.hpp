#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>
#include <string>

// Regression basis functions for LSMC.
// Under rough vol, the state is (S, v, Y, t) rather than just (S, t).
// We expose several basis choices so the caller can experiment.
class LSMCBasis {
public:
    enum Type {
        Monomial,     // 1, x, x^2, x^3, ... in each state variable
        Laguerre,     // weighted Laguerre polynomials L_0..L_k(x) — standard LS paper choice
        Chebyshev     // Chebyshev T_0..T_k(x) mapped to [-1,1] — better conditioning
    };

    // Build the feature vector for a single observation.
    // state = {S, v, Y} (normalised internally).
    // degree = polynomial degree per feature group.
    // includeVol / includeVolterra: whether to include v and Y as regressors.
    // For standard GBM: includeVol=false, includeVolterra=false → matches original LS paper.
    // For rough vol:    includeVol=true,  includeVolterra=true  → augmented state space.
    static std::vector<double> build(
        double S, double v, double Y,
        double S0, double xi0,   // normalisation constants
        int    degree,
        Type   type,
        bool   includeVol      = false,
        bool   includeVolterra = false)
    {
        std::vector<double> features;
        features.push_back(1.0);  // intercept always first

        // Normalise state variables to improve regression conditioning
        const double s  = S / S0;         // dimensionless spot
        const double sv = std::sqrt(std::max(v, 0.0) / std::max(xi0, 1e-12)); // vol/volATM
        const double sy = Y;              // Volterra is already dimensionless

        appendPolynomials(features, s,  degree, type);
        if (includeVol)
            appendPolynomials(features, sv, degree, type);
        if (includeVolterra)
            appendPolynomials(features, sy, degree, type);

        // Cross terms S*v and S*Y add meaningful interaction for rough vol
        if (includeVol && degree >= 2)
            features.push_back(s * sv);
        if (includeVolterra && degree >= 2)
            features.push_back(s * sy);

        return features;
    }

    // Batch version: builds design matrix X[nObs][nFeatures]
    static std::vector<std::vector<double>> buildMatrix(
        const std::vector<double>& S_vec,
        const std::vector<double>& v_vec,
        const std::vector<double>& Y_vec,
        double S0, double xi0,
        int degree, Type type,
        bool includeVol, bool includeVolterra)
    {
        const int n = static_cast<int>(S_vec.size());
        std::vector<std::vector<double>> X(n);
        for (int i = 0; i < n; ++i)
            X[i] = build(S_vec[i], v_vec[i], Y_vec[i],
                         S0, xi0, degree, type,
                         includeVol, includeVolterra);
        return X;
    }

    static int featureCount(int degree, Type type,
                            bool includeVol, bool includeVolterra)
    {
        int groups = 1;  // S always
        if (includeVol)      ++groups;
        if (includeVolterra) ++groups;
        int count = 1 + groups * degree;  // 1 intercept + degree terms per group
        if (includeVol      && degree >= 2) ++count;  // S*v cross
        if (includeVolterra && degree >= 2) ++count;  // S*Y cross
        return count;
    }

private:
    static void appendPolynomials(std::vector<double>& out,
                                  double x, int degree, Type type)
    {
        switch (type) {
            case Monomial:
                for (int k = 1; k <= degree; ++k)
                    out.push_back(std::pow(x, static_cast<double>(k)));
                break;

            case Laguerre:
                // Weighted Laguerre: e^{-x/2} * L_k(x) for x >= 0
                // Clamp negative x to 0 (can happen with rough vol)
                {
                    const double xp = std::max(x, 0.0);
                    const double w  = std::exp(-0.5 * xp);
                    out.push_back(w * laguerreL(1, xp));
                    for (int k = 2; k <= degree; ++k)
                        out.push_back(w * laguerreL(k, xp));
                }
                break;

            case Chebyshev:
                // Map x to [-1,1] via tanh for unbounded state variables
                {
                    const double z = std::tanh(x);
                    out.push_back(z);                // T_1(z) = z
                    if (degree >= 2) out.push_back(2.0*z*z - 1.0);      // T_2
                    if (degree >= 3) out.push_back(4.0*z*z*z - 3.0*z);  // T_3
                    for (int k = 4; k <= degree; ++k) {
                        // T_k = 2z*T_{k-1} - T_{k-2}
                        double tkm2 = out[out.size() - 2];
                        double tkm1 = out[out.size() - 1];
                        out.push_back(2.0 * z * tkm1 - tkm2);
                    }
                }
                break;
        }
    }

    // Laguerre polynomial L_k(x) via recurrence
    static double laguerreL(int k, double x) {
        if (k == 0) return 1.0;
        if (k == 1) return 1.0 - x;
        double lk2 = 1.0, lk1 = 1.0 - x;
        for (int j = 2; j <= k; ++j) {
            double lk = ((2*j - 1 - x) * lk1 - (j-1) * lk2) / j;
            lk2 = lk1;
            lk1 = lk;
        }
        return lk1;
    }
};