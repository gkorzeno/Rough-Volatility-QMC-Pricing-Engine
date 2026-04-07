// src/pricing/LongstaffSchwartz.hpp
#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include "../payoffs/Payoff.hpp"

class LongstaffSchwartz {
    // Fit a + b*x + c*x^2 to (x, y) pairs — basis for continuation value
    static std::vector<double> fitPolynomial(
        const std::vector<double>& x,
        const std::vector<double>& y)
    {
        // Simple least squares for degree-2 polynomial
        int n = x.size();
        double s0=n, s1=0, s2=0, s3=0, s4=0;
        double r0=0, r1=0, r2=0;

        for (int i = 0; i < n; i++) {
            double xi = x[i];
            s1 += xi; s2 += xi*xi; s3 += xi*xi*xi; s4 += xi*xi*xi*xi;
            r0 += y[i]; r1 += y[i]*xi; r2 += y[i]*xi*xi;
        }

        // Solve 3x3 normal equations (hardcoded for degree 2)
        // [s0 s1 s2] [a]   [r0]
        // [s1 s2 s3] [b] = [r1]
        // [s2 s3 s4] [c]   [r2]
        // Use Cramer's rule or simple Gaussian elimination
        // (for production use Eigen — this is illustrative)
        double A[3][4] = {
            {s0, s1, s2, r0},
            {s1, s2, s3, r1},
            {s2, s3, s4, r2}
        };

        // Gaussian elimination
        for (int col = 0; col < 3; col++) {
            for (int row = col+1; row < 3; row++) {
                double factor = A[row][col] / A[col][col];
                for (int k = col; k <= 3; k++)
                    A[row][k] -= factor * A[col][k];
            }
        }
        double c = A[2][3] / A[2][2];
        double b = (A[1][3] - A[1][2]*c) / A[1][1];
        double a = (A[0][3] - A[0][2]*c - A[0][1]*b) / A[0][0];

        return {a, b, c};
    }

public:
    static double price(
        const std::vector<std::vector<double>>& paths,
        const Payoff& payoff,
        double r,
        double dt)
    {
        int numPaths = paths.size();
        int steps    = paths[0].size() - 1;
        double disc  = std::exp(-r * dt);

        // Cashflow matrix — when does each path exercise?
        std::vector<double> cashflow(numPaths);
        for (int i = 0; i < numPaths; i++)
            cashflow[i] = payoff(paths[i][steps]);  // terminal payoff

        // Backward induction
        for (int t = steps - 1; t >= 1; t--) {
            // Find in-the-money paths at time t
            std::vector<int>    itmIdx;
            std::vector<double> S_itm, CV_itm;

            for (int i = 0; i < numPaths; i++) {
                double intrinsic = payoff(paths[i][t]);
                if (intrinsic > 0) {
                    itmIdx.push_back(i);
                    S_itm.push_back(paths[i][t]);
                    CV_itm.push_back(disc * cashflow[i]);  // discounted future cashflow
                }
            }

            if (itmIdx.empty()) continue;

            // Regress continuation value on current stock price
            auto coeffs = fitPolynomial(S_itm, CV_itm);

            // Exercise decision: intrinsic > estimated continuation?
            for (int idx : itmIdx) {
                double S          = paths[idx][t];
                double intrinsic  = payoff(S);
                double contVal    = coeffs[0] + coeffs[1]*S + coeffs[2]*S*S;

                if (intrinsic >= contVal)
                    cashflow[idx] = intrinsic;  // exercise early
                else
                    cashflow[idx] *= disc;       // hold — discount forward cashflow
            }

            // // Discount all other paths
            // for (int i = 0; i < numPaths; i++) {
            //     std::vector<bool> isITM(numPaths,false);
            //     for(int idx: itmIdx) isITM[idx]=true;
            //     // bool isItm = std::find(itmIdx.begin(), itmIdx.end(), i) != itmIdx.end();
            //     // if (!isItm) cashflow[i] *= disc;
            //     if(!isITM[i]) cashflow[i] *= disc;
            // }

            // Replace the "Discount all other paths" block with this:
            std::vector<bool> isITM(numPaths, false);
            for (int idx : itmIdx) isITM[idx] = true;

            for (int i = 0; i < numPaths; i++) {
                if (!isITM[i]) cashflow[i] *= disc;
            }
        }

        double sum = 0.0;
        for (double cf : cashflow) sum += cf;
        return sum / numPaths;
    }
};