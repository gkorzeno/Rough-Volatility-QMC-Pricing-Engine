#pragma once
#include <cmath>
#include <functional>
#include <stdexcept>
#include <vector>
#include "../core/Cholesky.hpp"
#include "../core/Random.hpp"

class FractionalBrownianMotion {
public:
    static double covariance(double t, double s, double hurst) {
        if (hurst <= 0.0 || hurst >= 1.0)
            throw std::runtime_error("FractionalBrownianMotion: H must be in (0,1)");
        return 0.5 * (std::pow(t, 2.0 * hurst)
                    + std::pow(s, 2.0 * hurst)
                    - std::pow(std::abs(t - s), 2.0 * hurst));
    }

    static std::vector<double> timeGrid(double T, int steps) {
        if (T <= 0.0 || steps <= 0)
            throw std::runtime_error("FractionalBrownianMotion: invalid horizon/grid");

        std::vector<double> times(steps);
        const double dt = T / steps;
        for (int i = 0; i < steps; ++i)
            times[i] = (i + 1) * dt;
        return times;
    }

    static std::vector<std::vector<double>> covarianceMatrix(
        const std::vector<double>& times,
        double hurst,
        double jitter = 1e-12)
    {
        const int n = static_cast<int>(times.size());
        std::vector<std::vector<double>> cov(n, std::vector<double>(n, 0.0));

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j <= i; ++j) {
                cov[i][j] = covariance(times[i], times[j], hurst);
                cov[j][i] = cov[i][j];
            }
            cov[i][i] += jitter;
        }
        return cov;
    }

    static std::vector<double> samplePath(
        double T,
        int steps,
        double hurst,
        Random& rng)
    {
        const std::vector<double> times = timeGrid(T, steps);
        const std::vector<std::vector<double>> cov = covarianceMatrix(times, hurst);
        const std::vector<std::vector<double>> L = Cholesky::decompose(cov);

        std::vector<double> z(steps, 0.0);
        for (int i = 0; i < steps; ++i)
            z[i] = rng.normal();

        std::vector<double> path(steps, 0.0);
        for (int i = 0; i < steps; ++i) {
            for (int j = 0; j <= i; ++j)
                path[i] += L[i][j] * z[j];
        }
        return path;
    }

    static std::vector<double> increments(const std::vector<double>& path) {
        std::vector<double> increments(path.size(), 0.0);
        double prev = 0.0;
        for (std::size_t i = 0; i < path.size(); ++i) {
            increments[i] = path[i] - prev;
            prev = path[i];
        }
        return increments;
    }
};
