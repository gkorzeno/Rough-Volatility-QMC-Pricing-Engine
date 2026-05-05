#pragma once
#include <vector>
#include <cmath>
#include "../core/Random.hpp"
#include "../core/Cholesky.hpp"
#include "../stochasticProcess/MultiDimensionalProcess.hpp"

class MultiDimensionalEulerMaruyama {
public:
    static std::vector<double> stepWithZ(
        const MultiDimensionalProcess& process,
        const std::vector<double>& x,
        double t,
        double dt,
        const std::vector<double>& z)
    {
        const int dim = static_cast<int>(x.size());
        if (static_cast<int>(z.size()) != dim)
            throw std::runtime_error("MultiDimensionalEulerMaruyama::stepWithZ dimension mismatch");

        auto mu = process.drift(x, t);
        auto sigma = process.diffusion(x, t);
        auto L = process.diffusionIsCholesky()
            ? sigma
            : Cholesky::decompose(sigma);

        std::vector<double> next(dim);
        const double sqrtDt = std::sqrt(dt);
        for (int i = 0; i < dim; i++) {
            double noise = 0.0;
            for (int j = 0; j <= i; j++)
                noise += L[i][j] * z[j];
            next[i] = x[i] + mu[i] * dt + noise * sqrtDt;
        }
        return next;
    }

    static std::vector<double> step(
        const MultiDimensionalProcess& process,
        const std::vector<double>& x,
        double t,
        double dt,
        Random& rng)
    {
        int dim = x.size();

        std::vector<double> z(dim);
        for (int i = 0; i < dim; i++) z[i] = rng.normal();
        return stepWithZ(process, x, t, dt, z);
    }
};
