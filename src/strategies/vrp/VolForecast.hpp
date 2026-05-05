#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>
#include "../../roughvol/RoughBergomiModel.hpp"
#include "../../surface/VolSurface.hpp"

// Estimates forward realized volatility using rBergomi Monte Carlo.
// Returns annualised vol as the sqrt of expected variance under the model.
class VolForecast {
public:
    struct Forecast {
        double meanRealizedVol;   // E[sqrt(1/T integral v_t dt)]
        double stdRealizedVol;    // cross-path std dev
        double volOfVol;          // std dev of variance (not vol)
        int    paths;
    };

    // Simulate rBergomi paths, compute realised variance on each,
    // return the distribution of annualised realised vol.
    static Forecast forecast(
        const RoughBergomiModel::Parameters& params,
        double spot,
        double rate,
        double horizon,           // years
        int    steps,
        int    paths,
        RoughBergomiModel::SamplingMethod method = RoughBergomiModel::RQMC,
        const std::string& dirFile = "docs/new-joe-kuo-6.21201.txt",
        std::uint64_t seed = 42)
    {
        if (paths <= 0 || steps <= 0)
            throw std::runtime_error("VolForecast: paths and steps must be positive");

        const double dt = horizon / steps;
        Random rng(static_cast<unsigned int>(seed));

        std::vector<double> realizedVols;
        realizedVols.reserve(paths);

        for (int i = 0; i < paths; ++i) {
            // Simulate a full path and collect variance at each step
            RoughBergomiModel::PathResult path =
                RoughBergomiModel::simulatePath(params, spot, rate, horizon, steps, rng);

            // Realized variance = time-average of instantaneous variance
            double sumVar = 0.0;
            for (int j = 1; j <= steps; ++j)
                sumVar += std::max(path.variance[j], 0.0);

            const double realizedVar = sumVar / steps;
            realizedVols.push_back(std::sqrt(realizedVar));
        }

        // Statistics across paths
        double mean = 0.0;
        for (double v : realizedVols) mean += v;
        mean /= paths;

        double sq = 0.0, sqVar = 0.0;
        for (double v : realizedVols) {
            sq    += (v - mean) * (v - mean);
            sqVar += (v * v - mean * mean) * (v * v - mean * mean);
        }

        return {mean,
                std::sqrt(sq / std::max(1, paths - 1)),
                std::sqrt(sqVar / std::max(1, paths - 1)),
                paths};
    }

    // Convenience: extract the ATM implied vol from a surface at a given expiry
    static double atmImpliedVol(const VolSurface& surface,
                                double spot, double expiry) {
        return surface.impliedVol(expiry, spot);
    }

    // VRP = implied vol - forecasted realized vol (positive = options expensive)
    static double vrp(double impliedVol, double forecastedVol) {
        return impliedVol - forecastedVol;
    }
};