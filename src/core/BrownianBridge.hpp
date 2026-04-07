// src/core/BrownianBridge.hpp
#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>

class BrownianBridge {
public:
    // Input:  Z[0..steps-1] — standard normals (e.g. from probit(sobol))
    // Output: dW[0..steps-1] — Brownian increments ~ N(0, dt) each
    //         constructed via Brownian bridge bisection so that
    //         the best Sobol dimensions drive the largest-scale moves
    static std::vector<double> build(
        const std::vector<double>& Z,
        double dt)
    {
        int steps = Z.size();
        if (steps == 0) throw std::runtime_error("BrownianBridge: empty input");

        // W[0] = 0 (fixed), W[steps] = terminal
        // We fill W[0..steps] then diff to get increments
        std::vector<double> W(steps + 1, 0.0);

        // --- Step 1: terminal value uses Z[0] (best Sobol dimension) ---
        W[steps] = std::sqrt(static_cast<double>(steps) * dt) * Z[0];

        // --- Step 2: bisection using a proper index mapping ---
        // We walk through Z[1], Z[2], ... in bisection order.
        // Each entry in the queue is (left_index, right_index).
        // We process level by level so z-index increments cleanly.

        struct Segment { int l, r; };
        std::vector<Segment> current, next_level;
        current.push_back({0, steps});

        int zi = 1;  // next Z index to consume

        while (!current.empty() && zi < steps) {
            next_level.clear();

            for (auto& seg : current) {
                if (zi >= steps) break;

                int l   = seg.l;
                int r   = seg.r;
                int mid = (l + r) / 2;

                if (mid == l || mid == r) continue;  // segment too small

                double tl = static_cast<double>(l)   * dt;
                double tm = static_cast<double>(mid)  * dt;
                double tr = static_cast<double>(r)    * dt;

                // Conditional mean and std of W(tm) | W(tl), W(tr)
                double condMean = W[l] + (tm - tl) / (tr - tl) * (W[r] - W[l]);
                double condStd  = std::sqrt((tm - tl) * (tr - tm) / (tr - tl));

                W[mid] = condMean + condStd * Z[zi++];

                // Push children for next level (only if they have a midpoint)
                if (mid - l >= 2) next_level.push_back({l,   mid});
                if (r - mid >= 2) next_level.push_back({mid, r  });
            }

            current = next_level;
        }

        // Any W[i] still 0 that wasn't reached by bisection: fill by linear interp
        // (only happens when steps is not a power of 2)
        for (int i = 1; i < steps; i++) {
            if (W[i] == 0.0 && i != steps) {
                // find nearest filled neighbours
                int l = i - 1;
                int r = i + 1;
                while (r <= steps && W[r] == 0.0 && r != steps) r++;
                // linear interpolation between W[l] and W[r]
                for (int k = l + 1; k < r; k++) {
                    double alpha = static_cast<double>(k - l) / (r - l);
                    W[k] = W[l] + alpha * (W[r] - W[l]);
                }
                i = r - 1;  // skip ahead
            }
        }

        // Convert Brownian path to increments
        std::vector<double> dW(steps);
        for (int i = 0; i < steps; i++)
            dW[i] = W[i + 1] - W[i];   // ~ N(0, dt) marginally

        return dW;
    }
};