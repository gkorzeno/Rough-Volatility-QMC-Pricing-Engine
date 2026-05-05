// src/control/LinearPolicy.hpp
#pragma once
#include "Policy.hpp"
#include <cmath>

// u(x,t) = -K * x  — standard LQR form
// For controlled OU, optimal K = (-theta + sqrt(theta^2 + sigma^2/r))
// where r is the control cost weight
class LinearPolicy : public Policy {
    double K;
public:
    LinearPolicy(double K_) : K(K_) {}
    double control(double x, double t) const override { return -K * x; }

    // Compute optimal K for OU with control cost r
    static double optimalGain(double theta, double sigma, double r) {
        return -theta + std::sqrt(theta*theta + sigma*sigma / r);
    }
};