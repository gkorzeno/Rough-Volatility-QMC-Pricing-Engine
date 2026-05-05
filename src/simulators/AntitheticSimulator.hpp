// // src/simulators/AntitheticSimulator.hpp
// #pragma once
// #include <cmath>
// #include <vector>
// #include "../integrators/EulerMaruyama.hpp"

// template<typename Integrator>
// class AntitheticSimulator {
// public:
//     static std::vector<double> simulate(
//         const StochasticProcess& process,
//         double x0, double T, double dt, int paths)
//     {
//         // paths must be even — each pair shares a noise sequence
//         int steps = static_cast<int>(T / dt);
//         std::vector<double> results(paths);
//         Random rng;

//         for (int i = 0; i < paths; i += 2) {
//             std::vector<double> zs(steps);
//             for (int j = 0; j < steps; j++) zs[j] = rng.normal();

//             // Original path
//             double x1 = x0, t = 0.0;
//             for (int j = 0; j < steps; j++) {
//                 double mu  = process.drift(x1, t);
//                 double sig = process.diffusion(x1, t);
//                 x1 += mu * dt + sig * std::sqrt(dt) * zs[j];
//                 t  += dt;
//             }

//             // Antithetic path — same noise negated
//             double x2 = x0; t = 0.0;
//             for (int j = 0; j < steps; j++) {
//                 double mu  = process.drift(x2, t);
//                 double sig = process.diffusion(x2, t);
//                 x2 += mu * dt + sig * std::sqrt(dt) * (-zs[j]);
//                 t  += dt;
//             }

//             results[i]   = x1;
//             results[i+1] = x2;
//         }
//         return results;
//     }
// };

#pragma once
#include <cmath>
#include <vector>
#include <stdexcept>
#include "../core/Random.hpp"
#include "../stochasticProcess/StochasticProcess.hpp"

template<typename Integrator>
class AntitheticSimulator {
public:
    static std::vector<double> simulate(
        const StochasticProcess& process,
        double x0, double T, double dt, int paths)

    {

        if (paths % 2 != 0) {
            throw std::runtime_error("AntitheticSimulator requires an even number of paths");
        }

        int steps = static_cast<int>(T / dt);
        std::vector<double> results(paths);
        Random rng;

        for (int i = 0; i < paths; i += 2) {

            double x1 = x0;
            double x2 = x0;
            double t  = 0.0;

            for (int j = 0; j < steps; j++) {

                double z = rng.normal();

                x1 = Integrator::stepWithZ(process, x1, t, dt,  z);
                x2 = Integrator::stepWithZ(process, x2, t, dt, -z);

                t += dt;
            }

            results[i]   = x1;
            results[i+1] = x2;
            // results[i/2] = 0.5 * (x1 + x2);
        }

        return results;
    }
};