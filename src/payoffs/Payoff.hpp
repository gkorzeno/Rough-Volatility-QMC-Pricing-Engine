// src/payoffs/Payoff.hpp
#pragma once

class Payoff {
public:
    virtual double operator()(double S) const = 0;
    virtual ~Payoff() {}
};