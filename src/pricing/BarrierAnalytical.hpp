// src/pricing/BarrierAnalytical.hpp
#pragma once
#include <cmath>
#include <stdexcept>
#include "../barriers/Barrier.hpp"

class BarrierAnalytical {
    static double N(double x) {
        return 0.5 * std::erfc(-x / std::sqrt(2.0));
    }

    // Core terms used across all barrier formulas
    struct Terms {
        double phi, eta;   // direction parameters
        double d1, d2, d3, d4;
        double A, B, C, D;
    };

public:
    // Exact BS price for European barrier options
    // Reference: Haug "Complete Guide to Option Pricing Formulas" Ch.4
    static double price(
        double S, double K, double H,
        double r, double sigma, double T,
        BarrierType type)
    {
        double mu    = (r - 0.5*sigma*sigma) / (sigma*sigma);
        double lam   = std::sqrt(mu*mu + 2*r/(sigma*sigma));
        double sqrtT = std::sqrt(T);

        double x1 = std::log(S/K)  / (sigma*sqrtT) + (1+mu)*sigma*sqrtT;
        double x2 = std::log(S/H)  / (sigma*sqrtT) + (1+mu)*sigma*sqrtT;
        double y1 = std::log(H*H/(S*K)) / (sigma*sqrtT) + (1+mu)*sigma*sqrtT;
        double y2 = std::log(H/S) / (sigma*sqrtT) + (1+mu)*sigma*sqrtT;

        double A  =  S * N(x1)
                   - K * std::exp(-r*T) * N(x1 - sigma*sqrtT);
        double B  =  S * N(x2)
                   - K * std::exp(-r*T) * N(x2 - sigma*sqrtT);
        double C  =  S * std::pow(H/S, 2*(mu+1)) * N(y1)
                   - K * std::exp(-r*T) * std::pow(H/S, 2*mu) * N(y1 - sigma*sqrtT);
        double D  =  S * std::pow(H/S, 2*(mu+1)) * N(y2)
                   - K * std::exp(-r*T) * std::pow(H/S, 2*mu) * N(y2 - sigma*sqrtT);

        switch (type) {
            case BarrierType::DownAndOut: return (S >= H) ? A - C : 0.0;
            // K <= H branch for up barrier calls: Reiner-Rubinstein / Haug form.
            // The prior code had a sign error on the image terms (C, D).
            case BarrierType::UpAndOut:   return (S <= H) ? A - B - C + D : 0.0;
            case BarrierType::DownAndIn:  return (S >= H) ? C : A;
            case BarrierType::UpAndIn:    return (S <= H) ? B + C - D : A;
            default: throw std::runtime_error("Unknown barrier type");
        }
    }
};