// src/pricing/BlackScholes.hpp
#pragma once
#include "../greeks/Greeks.hpp"
#include <cmath>


class BlackScholes {
    static double N(double x) {
        return 0.5 * std::erfc(-x / std::sqrt(2.0));
    }
public:
    static double callPrice(double S, double K, double r,
                            double sigma, double T) {
        double d1 = (std::log(S/K) + (r + 0.5*sigma*sigma)*T)
                    / (sigma * std::sqrt(T));
        double d2 = d1 - sigma * std::sqrt(T);
        return S * N(d1) - K * std::exp(-r*T) * N(d2);
    }

    static double putPrice(double S, double K, double r,
                           double sigma, double T) {
        return callPrice(S, K, r, sigma, T) - S + K * std::exp(-r * T);
    }

    static Greeks greeks(
        double S, double K, double r, double sigma, double T)
    {
        double d1 = (std::log(S/K) + (r + 0.5*sigma*sigma)*T)
                / (sigma * std::sqrt(T));
        double d2 = d1 - sigma * std::sqrt(T);
        double disc = std::exp(-r * T);
        double nd1  = N(d1);
        double nd2  = N(d2);
        double npd1 = std::exp(-0.5*d1*d1) / std::sqrt(2.0 * 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067);

        Greeks g;
        g.price = S * nd1 - K * disc * nd2;
        g.delta = nd1;
        g.gamma = npd1 / (S * sigma * std::sqrt(T));
        g.vega  = S * npd1 * std::sqrt(T);
        g.theta = (-S * npd1 * sigma / (2*std::sqrt(T))
               - r * K * disc * nd2) / 365.0;
        g.rho   = K * T * disc * nd2;
        return g;
    }

    static double vega(double S, double K, double r, double sigma, double T) {
        double d1  = (std::log(S/K) + (r + 0.5*sigma*sigma)*T)
                 / (sigma * std::sqrt(T));
        double npd1 = std::exp(-0.5*d1*d1) / std::sqrt(2.0 * 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067);
        return S * npd1 * std::sqrt(T);
    }
};
