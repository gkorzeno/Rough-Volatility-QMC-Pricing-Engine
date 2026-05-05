// src/payoffs/MultiAssetPayoff.hpp
#pragma once
#include <vector>
#include <algorithm>
#include <numeric>
#include <stdexcept>

class MultiAssetPayoff {
public:
    virtual double operator()(const std::vector<double>& S) const = 0;
    virtual ~MultiAssetPayoff() {}
};

// Max(w1*S1 + w2*S2 + ... - K, 0)
class BasketCallPayoff : public MultiAssetPayoff {
    std::vector<double> weights;
    double K;
public:
    BasketCallPayoff(const std::vector<double>& weights_, double K_)
        : weights(weights_), K(K_) {}

    double operator()(const std::vector<double>& S) const override {
        double basket = 0.0;
        for (size_t i = 0; i < S.size(); i++)
            basket += weights[i] * S[i];
        return std::max(basket - K, 0.0);
    }
};

// Max(S1 - S2 - K, 0) — Margrabe spread option
class SpreadCallPayoff : public MultiAssetPayoff {
    double K;
    int i, j;   // which two assets to spread
public:
    SpreadCallPayoff(double K_, int i_ = 0, int j_ = 1)
        : K(K_), i(i_), j(j_) {}

    double operator()(const std::vector<double>& S) const override {
        return std::max(S[i] - S[j] - K, 0.0);
    }
};

// Max(max(S1,...,Sn) - K, 0)
class RainbowCallPayoff : public MultiAssetPayoff {
    double K;
public:
    RainbowCallPayoff(double K_) : K(K_) {}

    double operator()(const std::vector<double>& S) const override {
        double best = *std::max_element(S.begin(), S.end());
        return std::max(best - K, 0.0);
    }
};

// Max(K - min(S1,...,Sn), 0) — worst-of put
class WorstOfPutPayoff : public MultiAssetPayoff {
    double K;
public:
    WorstOfPutPayoff(double K_) : K(K_) {}

    double operator()(const std::vector<double>& S) const override {
        double worst = *std::min_element(S.begin(), S.end());
        return std::max(K - worst, 0.0);
    }
};