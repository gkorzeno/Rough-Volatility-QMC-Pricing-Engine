#pragma once
#include <random>

class Random {
private:
    std::mt19937 generator;
    std::normal_distribution<double> normal_dist;
    std::uniform_real_distribution<double> uniform_dist;

public:
    Random()
        : generator(std::random_device{}()),
          normal_dist(0.0, 1.0),
          uniform_dist(0.0, 1.0) {}

    Random(unsigned int seed) : generator(seed),
        normal_dist(0.0, 1.0),
        uniform_dist(0.0, 1.0) {}

    double normal() {
        return normal_dist(generator);
        // return std::uniform_real_distribution<double>(0.0, 1.0)(generator);
    }

    double uniform() {
        return uniform_dist(generator);
    }

    int poisson(double lambda) {
        // Knuth algorithm — efficient for small lambda * dt (which is typical)
        //Knuth algorithm could slow down high-frequency jumps
        //Potential Upgrade: 
        // std::poisson_distribution<int> poisson_dist(lambda);
        // return poisson_dist(generator);
        double L = std::exp(-lambda);
        int k = 0;
        double p = 1.0;
        do {
            k++;
            p *= uniform(); // your existing U(0,1) sampler
        } while (p > L);
        return k - 1;
    }
};