// src/stochasticProcess/CorrelatedGBM.hpp
#pragma once
#include "MultiDimensionalProcess.hpp"
#include "../core/Cholesky.hpp"
#include "CholeskyProcess.hpp"
#include <vector>
#include <cmath>
#include <stdexcept>

// Models N correlated GBMs
// State vector x[i] = S_i (asset prices)
class CorrelatedGBM : public CholeskyProcess {
private:
    std::vector<double> mu;        // drift per asset
    std::vector<double> sigma;     // volatility per asset
    std::vector<std::vector<double>> correlationMatrix;
    std::vector<std::vector<double>> L;  // Cholesky factor of correlation matrix
    int n;

public:
    CorrelatedGBM(
        const std::vector<double>& mu_,
        const std::vector<double>& sigma_,
        const std::vector<std::vector<double>>& correlationMatrix)
        : mu(mu_), sigma(sigma_), correlationMatrix(correlationMatrix), n(mu_.size())
    {
        if (mu_.size() != sigma_.size() ||
            correlationMatrix.size() != (size_t)n)
            throw std::runtime_error("CorrelatedGBM: dimension mismatch");

        L = Cholesky::decompose(correlationMatrix);
    }

    std::vector<double> drift(
        const std::vector<double>& x, double t) const override
    {
        std::vector<double> d(n);
        for (int i = 0; i < n; i++)
            d[i] = mu[i] * x[i];
        return d;
    }

    // Returns pre-factored L scaled by sigma_i * S_i
    // Row i: L[i][j] * sigma[i] * x[i]
    std::vector<std::vector<double>> diffusion(
        const std::vector<double>& x, double t) const override
    {
        std::vector<std::vector<double>> D(n, std::vector<double>(n, 0.0));
        for (int i = 0; i < n; i++)
            for (int j = 0; j <= i; j++)
                D[i][j] = L[i][j] * sigma[i] * x[i];
        return D;
    }

    int size() const { return n; }
    const std::vector<double>& drifts() const { return mu; }
    const std::vector<double>& volatilities() const { return sigma; }
    const std::vector<std::vector<double>>& correlation() const { return correlationMatrix; }
    const std::vector<std::vector<double>>& choleskyFactor() const { return L; }

    // Risk-neutral drift override (for pricing)
    void setRiskNeutral(double r) {
        for (int i = 0; i < n; i++)
            mu[i] = r;
    }
};
