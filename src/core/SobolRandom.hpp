// src/core/SobolRandom.hpp
#pragma once
#include <cmath>
#include <vector>
#include <numbers>    // C++20 — has std::numbers::pi
#include "SobolSequence.hpp"

class SobolRandom {
private:
    static constexpr double PI = 3.14159265358979323846;

    SobolSequence sobol;
    std::vector<double> buffer;
    int bufferIdx;

    std::pair<double,double> boxMuller(double u1, double u2) {
        double r     = std::sqrt(-2.0 * std::log(u1 + 1e-300));
        double theta = 2.0 * PI * u2;
        return {r * std::cos(theta), r * std::sin(theta)};
    }

public:
    // SobolRandom() : sobol(2), bufferIdx(2) {}
    SobolRandom()
        : sobol(2, "docs/new-joe-kuo-6.21201.txt"), bufferIdx(2) {}

    double normal() {
        if (bufferIdx >= 2) {
            auto uv       = sobol.next();
            auto [z1, z2] = boxMuller(uv[0], uv[1]);
            buffer        = {z1, z2};
            bufferIdx     = 0;
        }
        return buffer[bufferIdx++];
    }
};