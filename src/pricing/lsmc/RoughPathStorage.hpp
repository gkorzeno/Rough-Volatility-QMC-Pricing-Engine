#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>
#include "../../roughvol/RoughBergomiModel.hpp"
#include "../../core/Random.hpp"

// Stores full rBergomi simulation paths needed for LSMC.
// Under rough Bergomi the state is non-Markovian: the continuation
// value at time t_j depends on (S_j, v_j, Y_j) where Y_j is the
// Volterra process integral. We store all three.
struct RoughPathData {
    // [path][step]: spot price
    std::vector<std::vector<double>> S;
    // [path][step]: instantaneous variance v_t = xi0 * exp(eta*Y - 0.5*eta^2*t^{2H})
    std::vector<std::vector<double>> v;
    // [path][step]: Volterra process Y_t (the rough driver)
    std::vector<std::vector<double>> Y;
    // time grid t[0..nSteps]
    std::vector<double> times;

    int nPaths = 0;
    int nSteps = 0;
    double maturity = 0.0;
};

class RoughPathStorage {
public:
    // Simulate nPaths full rBergomi paths on a uniform grid of nSteps steps.
    // Uses the exact Volterra simulation from RoughBergomiModel::simulatePath.
    static RoughPathData simulate(
        const RoughBergomiModel::Parameters& params,
        double spot,
        double rate,
        double maturity,
        int    nSteps,
        int    nPaths,
        unsigned int seed = 42)
    {
        if (nSteps <= 0 || nPaths <= 0)
            throw std::runtime_error("RoughPathStorage: nSteps and nPaths must be positive");

        RoughPathData data;
        data.nPaths   = nPaths;
        data.nSteps   = nSteps;
        data.maturity = maturity;

        const double dt = maturity / nSteps;
        data.times.resize(nSteps + 1);
        for (int j = 0; j <= nSteps; ++j)
            data.times[j] = j * dt;

        data.S.assign(nPaths, std::vector<double>(nSteps + 1, 0.0));
        data.v.assign(nPaths, std::vector<double>(nSteps + 1, params.xi0));
        data.Y.assign(nPaths, std::vector<double>(nSteps + 1, 0.0));

        Random rng(seed);

        for (int i = 0; i < nPaths; ++i) {
            auto path = RoughBergomiModel::simulatePath(
                params, spot, rate, maturity, nSteps, rng);

            for (int j = 0; j <= nSteps; ++j) {
                data.S[i][j] = path.spot[j];
                data.v[i][j] = path.variance[j];
                data.Y[i][j] = path.volterra[j];
            }
        }

        return data;
    }

    // Extract a cross-section at step j: returns (S_j, v_j, Y_j) for all paths
    static void crossSection(
        const RoughPathData& data,
        int step,
        std::vector<double>& S_out,
        std::vector<double>& v_out,
        std::vector<double>& Y_out)
    {
        const int n = data.nPaths;
        S_out.resize(n);
        v_out.resize(n);
        Y_out.resize(n);
        for (int i = 0; i < n; ++i) {
            S_out[i] = data.S[i][step];
            v_out[i] = data.v[i][step];
            Y_out[i] = data.Y[i][step];
        }
    }
};