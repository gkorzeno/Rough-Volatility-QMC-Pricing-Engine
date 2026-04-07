// tests/test_sobol.cpp
#include <iostream>
#include <cassert>
#include <cmath>
#include "../src/core/SobolSequence.hpp"

const std::string DIR_FILE = "docs/new-joe-kuo-6.21201.txt";

void test_mean() {
    SobolSequence sobol(1, DIR_FILE);
    int N = 10000;
    double sum = 0.0;
    for (int i = 0; i < N; i++)
        sum += sobol.next()[0];

    double mean = sum / N;
    std::cout << "[Mean Test]     Mean = " << mean
              << "  (expected ~0.5)\n";
    assert(std::abs(mean - 0.5) < 1e-4 && "Mean should be ~0.5");
}

void test_variance() {
    SobolSequence sobol(1, DIR_FILE);
    int N = 10000;

    // Must skip the same N points first to get a fresh sequence
    // OR reset — depends on your API. Here we just construct fresh.
    double sum  = 0.0;
    double sum2 = 0.0;
    for (int i = 0; i < N; i++) {
        double x = sobol.next()[0];
        sum  += x;
        sum2 += x * x;
    }
    double mean = sum / N;
    double var  = sum2 / N - mean * mean;   // E[X^2] - E[X]^2

    std::cout << "[Variance Test] Variance = " << var
              << "  (expected ~0.0833 = 1/12)\n";
    assert(std::abs(var - 1.0/12.0) < 1e-4 && "Variance should be ~1/12");
}

void test_2d_uniformity() {
    SobolSequence sobol(2, DIR_FILE);
    std::cout << "[2D Uniformity] First 10 points:\n";

    // Known correct Sobol (0,2) sequence first 10 points
    // Dim 0: 0.5, 0.75, 0.25, 0.375, 0.875, 0.625, 0.125, ...
    // Dim 1: 0.5, 0.25, 0.75, 0.375, 0.875, 0.125, 0.625, ...
    for (int i = 0; i < 10; i++) {
        auto p = sobol.next();
        std::cout << "  " << p[0] << ", " << p[1] << "\n";
    }

    // Quantitative check: in [0,1]^2, all 4 quadrants should be hit
    // within the first 4 points of a proper Sobol sequence
    SobolSequence sobol2(2, DIR_FILE);
    bool q1=false, q2=false, q3=false, q4=false;
    for (int i = 0; i < 4; i++) {
        auto p = sobol2.next();
        bool hi_x = p[0] > 0.5;
        bool hi_y = p[1] > 0.5;
        if ( hi_x &&  hi_y) q1 = true;
        if (!hi_x &&  hi_y) q2 = true;
        if (!hi_x && !hi_y) q3 = true;
        if ( hi_x && !hi_y) q4 = true;
    }
    assert(q1 && q2 && q3 && q4 && "First 4 points must cover all 4 quadrants");
    std::cout << "  All 4 quadrants covered in first 4 points: PASS\n";
}

int main() {
    std::cout << "=== Sobol Sequence Tests ===\n\n";
    test_mean();
    test_variance();
    test_2d_uniformity();
    std::cout << "\nAll tests passed.\n";
    return 0;
}
