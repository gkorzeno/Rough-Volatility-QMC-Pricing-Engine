#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>
#include <numeric>

// Ordinary least squares regression for LSMC continuation value estimation.
// Solves beta = (X'X)^{-1} X'y via Cholesky decomposition.
// Keeps the regression numerically stable with Tikhonov regularisation.
class LSMCRegressor {
public:
    struct FitResult {
        std::vector<double> coefficients;
        double              rSquared;
        int                 nObs;
        int                 nFeatures;
    };

    // Fit: given design matrix X[n][p] and response y[n], solve for beta[p]
    static FitResult fit(
        const std::vector<std::vector<double>>& X,
        const std::vector<double>&              y,
        double lambda = 1e-6)  // Tikhonov regularisation
    {
        const int n = static_cast<int>(X.size());
        const int p = n > 0 ? static_cast<int>(X[0].size()) : 0;
        if (n < p)
            throw std::runtime_error("LSMCRegressor: underdetermined system (n < p)");

        // Build XtX [p x p] and Xty [p]
        std::vector<std::vector<double>> XtX(p, std::vector<double>(p, 0.0));
        std::vector<double>              Xty(p, 0.0);

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < p; ++j) {
                Xty[j] += X[i][j] * y[i];
                for (int k = 0; k <= j; ++k)
                    XtX[j][k] += X[i][j] * X[i][k];
            }
        }
        // Symmetrise and add regularisation on diagonal
        for (int j = 0; j < p; ++j) {
            XtX[j][j] += lambda;
            for (int k = 0; k < j; ++k)
                XtX[k][j] = XtX[j][k];
        }

        // Cholesky decomposition XtX = L L'
        auto L = choleskyDecompose(XtX);

        // Solve L L' beta = Xty via forward/back substitution
        auto beta = choleskysolve(L, Xty);

        // Compute R^2
        double yMean = std::accumulate(y.begin(), y.end(), 0.0) / n;
        double ssTot = 0.0, ssRes = 0.0;
        for (int i = 0; i < n; ++i) {
            double yHat = 0.0;
            for (int j = 0; j < p; ++j) yHat += X[i][j] * beta[j];
            ssRes += (y[i] - yHat) * (y[i] - yHat);
            ssTot += (y[i] - yMean) * (y[i] - yMean);
        }
        double r2 = (ssTot > 1e-12) ? 1.0 - ssRes / ssTot : 0.0;

        return {beta, r2, n, p};
    }

    // Predict: given fitted coefficients and a feature vector, return E[V|state]
    static double predict(const std::vector<double>& beta,
                          const std::vector<double>& features)
    {
        double val = 0.0;
        const int p = static_cast<int>(std::min(beta.size(), features.size()));
        for (int j = 0; j < p; ++j)
            val += beta[j] * features[j];
        return val;
    }

private:
    static std::vector<std::vector<double>> choleskyDecompose(
        const std::vector<std::vector<double>>& A)
    {
        const int n = static_cast<int>(A.size());
        std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j <= i; ++j) {
                double sum = 0.0;
                for (int k = 0; k < j; ++k)
                    sum += L[i][k] * L[j][k];
                if (i == j) {
                    double val = A[i][i] - sum;
                    L[i][j] = (val > 0.0) ? std::sqrt(val) : std::sqrt(1e-12);
                } else {
                    L[i][j] = (L[j][j] > 1e-12)
                        ? (A[i][j] - sum) / L[j][j]
                        : 0.0;
                }
            }
        }
        return L;
    }

    static std::vector<double> choleskysolve(
        const std::vector<std::vector<double>>& L,
        const std::vector<double>& b)
    {
        const int n = static_cast<int>(b.size());
        std::vector<double> z(n, 0.0), x(n, 0.0);

        // Forward substitution: L z = b
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < i; ++j) sum += L[i][j] * z[j];
            z[i] = (L[i][i] > 1e-12) ? (b[i] - sum) / L[i][i] : 0.0;
        }

        // Back substitution: L' x = z
        for (int i = n - 1; i >= 0; --i) {
            double sum = 0.0;
            for (int j = i + 1; j < n; ++j) sum += L[j][i] * x[j];
            x[i] = (L[i][i] > 1e-12) ? (z[i] - sum) / L[i][i] : 0.0;
        }

        return x;
    }
};