// src/payoffs/VanillaPayoffs.hpp
#pragma once
#include "Payoff.hpp"
#include <cmath>
#include <algorithm>

class CallPayoff : public Payoff {
    double K;
public:
    CallPayoff(double K_) : K(K_) {}
    double operator()(double S) const override { return std::max(S - K, 0.0); }
};

class PutPayoff : public Payoff {
    double K;
public:
    PutPayoff(double K_) : K(K_) {}
    double operator()(double S) const override { return std::max(K - S, 0.0); }
};

class DigitalCallPayoff : public Payoff {
    double K;
public:
    DigitalCallPayoff(double K_) : K(K_) {}
    double operator()(double S) const override { return S > K ? 1.0 : 0.0; }
};