// src/utils/Benchmark.hpp
#pragma once
#include <chrono>
#include <string>
#include <iostream>
#include <functional>

class Benchmark {
public:
    static void run(
        const std::string& label,
        std::function<void()> fn)
    {
        auto start = std::chrono::high_resolution_clock::now();
        fn();
        auto end   = std::chrono::high_resolution_clock::now();

        double ms = std::chrono::duration<double, std::milli>(end - start).count();
        std::cout << label << ": " << ms << " ms\n";
    }
};