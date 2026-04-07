// src/barriers/BarrierCrossing.hpp
#pragma once
#include <cmath>
#include <algorithm>
#include "../barriers/Barrier.hpp"

class BarrierCrossing {
public:
    // Probability that a GBM Brownian bridge crossed H between t and t+dt
    // given S(t) = s0 and S(t+dt) = s1
    // Formula: exp(-2 * log(s0/H) * log(s1/H) / (sigma^2 * dt))
    static double crossingProbability(
        double s0, double s1,
        double H,
        double sigma, double dt)
    {
        // double a = std::log(s0 / H);
        // double b = std::log(s1 / H);

        double x0 = std::log(s0);
        double x1 = std::log(s1);
        double h  = std::log(H);

        double a = x0 - h;
        double b = x1 - h;

        // If already on opposite sides, crossing is certain
        if (a * b <= 0.0) return 1.0;

        // If both on same side, use bridge probability
        double prob = std::exp(-2.0 * a * b / (sigma * sigma * dt));
        return std::min(prob, 1.0);
    }

    // Determine if barrier was crossed this step, using bridge detection
    // Returns true if a crossing occurred
    static bool crossedWithBridge(
        double s0, double s1,
        const Barrier& barrier,
        double sigma, double dt,
        double uniform)   // uniform random in [0,1] for probabilistic test
    {
        // double H = barrier.H;

        // // Discrete check: did we cross at observed points?
        // bool discreteCross = barrier.isUp()
        //     ? (s0 < H && s1 >= H) || (s0 >= H && s1 < H)
        //     : (s0 > H && s1 <= H) || (s0 <= H && s1 > H);

        // if (discreteCross) return true;

        // // Bridge check: could we have crossed between steps?
        // // Only check if both points are on the same side
        // bool s0above = (s0 >= H);
        // bool s1above = (s1 >= H);
        // if (s0above != s1above) return true;  // already caught above

        // // Both on same side — use bridge probability
        // double prob = crossingProbability(s0, s1, H, sigma, dt);
        // return uniform < prob;

        double H = barrier.H;

        bool s0above = (s0 >= H);
        bool s1above = (s1 >= H);

        // Definite cross — endpoints on opposite sides
        if (s0above != s1above) return true;

        // Both on same side — use bridge probability
        // Only apply bridge if the barrier is on the correct side
        // For up barrier: only matters if both points are BELOW H
        // For down barrier: only matters if both points are ABOVE H
        if (barrier.isUp() && s0above)   return false;  // both above, can't cross up
        if (!barrier.isUp() && !s0above) return false;  // both below, can't cross down

        double prob = crossingProbability(s0, s1, H, sigma, dt);
        return uniform < prob;
    }
};