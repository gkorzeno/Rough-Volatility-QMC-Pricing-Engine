// src/control/ConstantPolicy.hpp
#pragma once
#include "Policy.hpp"

class ConstantPolicy : public Policy {
    double u;
public:
    ConstantPolicy(double u_) : u(u_) {}
    double control(double x, double t) const override { return u; }
};