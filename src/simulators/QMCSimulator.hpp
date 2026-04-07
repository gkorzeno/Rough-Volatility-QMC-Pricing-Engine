// src/simulators/QMCSimulator.hpp
#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>
#include "../stochasticProcess/StochasticProcess.hpp"
#include "../core/SobolSequence.hpp"
#include "../core/BrownianBridge.hpp"

// Rational approximation coefficients for probit (Beasley-Springer-Moro)
static double probit(double u) {
    // Clamp to avoid log(0) or log(-inf)
    if (u <= 0.0) u = 1e-15;
    if (u >= 1.0) u = 1.0 - 1e-15;

    static const double a[] = {
        -3.969683028665376e+01,  2.209460984245205e+02,
        -2.759285104469687e+02,  1.383577518672690e+02,
        -3.066479806614716e+01,  2.506628277459239e+00
    };
    static const double b[] = {
        -5.447609879822406e+01,  1.615858368580409e+02,
        -1.556989798598866e+02,  6.680131188771972e+01,
        -1.328068155288572e+01
    };
    static const double c[] = {
        -7.784894002430293e-03, -3.223964580411365e-01,
        -2.400758277161838e+00, -2.549732539343734e+00,
         4.374664141464968e+00,  2.938163982698783e+00
    };
    static const double d[] = {
         7.784695709041462e-03,  3.224671290700398e-01,
         2.445134137142996e+00,  3.754408661907416e+00
    };

    const double pLow  = 0.02425;
    const double pHigh = 1.0 - pLow;

    double q, r;
    if (u < pLow) {
        q = std::sqrt(-2.0 * std::log(u));
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
               ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
    } else if (u <= pHigh) {
        q = u - 0.5;
        r = q * q;
        return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
               (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1.0);
    } else {
        q = std::sqrt(-2.0 * std::log(1.0 - u));
        return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
                ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
    }
}

template<typename Integrator>
class QMCSimulator {
public:
    static std::vector<double> simulate(
        const StochasticProcess& process,
        double x0, double T, double dt, int paths,
        const std::string& directionFile)
    {
        int steps = static_cast<int>(T / dt);
        std::vector<double> results(paths);

        // One Sobol dimension per timestep
        SobolSequence sobol(steps, directionFile);

        for (int i = 0; i < paths; i++) {
            // Get one quasi-random point in [0,1]^steps
            auto u = sobol.next();

            // Convert each uniform to a standard normal via probit (inverse CDF)
            // This avoids Box-Muller pairing issues entirely
            std::vector<double> z(steps);
            for (int j = 0; j < steps; j++)
                z[j] = probit(u[j]);

            // BrownianBridge reorders z-scores for better low-discrepancy coverage
            // It returns z-scores (not scaled increments), so pass z directly
            auto zBridge = BrownianBridge::build(z, dt);

            // Simulate path — stepWithZ expects a standard normal z, not dW
            double x = x0, t = 0.0;
            for (int j = 0; j < steps; j++) {
                // Extract z-score: bridge returns dW = z * sqrt(dt), so recover z
                double zj = zBridge[j] / std::sqrt(dt);
                x = Integrator::stepWithZ(process, x, t, dt, zj);
                t += dt;
            }
            results[i] = x;
        }
        return results;
    }
};