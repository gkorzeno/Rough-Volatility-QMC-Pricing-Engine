#pragma once
#include <vector>
#include "../core/Random.hpp"

class MultiDimensionalProcess {
public:

    // drift vector μ(x,t)
    virtual std::vector<double> drift(
        const std::vector<double>& x,
        double t) const = 0;

    // diffusion matrix σ(x,t)
    virtual std::vector<std::vector<double>> diffusion(
        const std::vector<double>& x,
        double t) const = 0;

    virtual std::vector<double> theoreticalMean(
        const std::vector<double>& x0, double t) const
    {
        return std::vector<double>(x0.size(), 0.0); // override per process
    }

    virtual std::vector<double> theoreticalVariance(double t) const {
        return std::vector<double>(); // override per process
    }

    virtual bool diffusionIsCholesky() const { return false; }

    virtual ~MultiDimensionalProcess() {}
};