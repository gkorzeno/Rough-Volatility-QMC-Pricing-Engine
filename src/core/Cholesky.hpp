// src/core/Cholesky.hpp
#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>

class Cholesky {
public:
    // Returns lower triangular L such that L*L^T = A
    static std::vector<std::vector<double>> decompose(
        const std::vector<std::vector<double>>& A)
    {
        int n = A.size();
        std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));

        for (int i = 0; i < n; i++) {
            for (int j = 0; j <= i; j++) {
                double sum = 0.0;
                for (int k = 0; k < j; k++)
                    sum += L[i][k] * L[j][k];

                if (i == j) {
                    double val = A[i][i] - sum;
                    if (val <= 0.0)
                        throw std::runtime_error("Matrix not positive definite");
                    L[i][j] = std::sqrt(val);
                } else {
                    L[i][j] = (A[i][j] - sum) / L[j][j];
                }
            }
        }
        return L;
    }
};