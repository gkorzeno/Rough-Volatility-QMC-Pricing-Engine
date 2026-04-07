// src/greeks/Greeks.hpp
#pragma once
#include <cstdio>

struct Greeks {
    double price;
    double delta;    // dV/dS
    double gamma;    // d²V/dS²
    double vega;     // dV/dsigma
    double theta;    // dV/dT  (per day)
    double rho;      // dV/dr

    void print() const {
        printf("Price: %8.4f\n", price);
        printf("Delta: %8.4f\n", delta);
        printf("Gamma: %8.4f\n", gamma);
        printf("Vega:  %8.4f\n", vega);
        printf("Theta: %8.4f\n", theta);
        printf("Rho:   %8.4f\n", rho);
    }
};