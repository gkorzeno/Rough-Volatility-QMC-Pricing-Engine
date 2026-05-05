// src/payoffs/BarrierPayoff.hpp
#pragma once
#include <cmath>
#include <algorithm>
#include "../payoffs/Payoff.hpp"
#include "../barriers/Barrier.hpp"

class BarrierCallPayoff : public Payoff {
    double K;
public:
    BarrierCallPayoff(double K_) : K(K_) {}
    double operator()(double S) const override {
        return std::max(S - K, 0.0);
    }
};

class BarrierPutPayoff : public Payoff {
    double K;
public:
    BarrierPutPayoff(double K_) : K(K_) {}
    double operator()(double S) const override {
        return std::max(K - S, 0.0);
    }
};