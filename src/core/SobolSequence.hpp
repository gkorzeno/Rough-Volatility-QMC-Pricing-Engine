// // src/core/SobolSequence.hpp
// #pragma once
// #include <vector>
// #include <cmath>
// #include <stdexcept>

// class SobolSequence {
// private:
//     int dim;
//     int index;
//     std::vector<std::vector<unsigned int>> directionNumbers;
//     std::vector<unsigned int> state;

//     static const int BITS = 32;

//     void initDirectionNumbers() {
//         // Direction numbers for first 2 dimensions (Joe & Kuo 2010 table)
//         // For production: load full table from file for high dimensions
//         directionNumbers.resize(dim);

//         // Dimension 1: all 1s shifted
//         directionNumbers[0].resize(BITS);
//         for (int i = 0; i < BITS; i++)
//             directionNumbers[0][i] = 1u << (BITS - 1 - i);

//         if (dim < 2) return;

//         // Dimension 2: s=1, a=0, m={1}
//         directionNumbers[1].resize(BITS);
//         directionNumbers[1][0] = 1u << (BITS - 1);
//         for (int i = 1; i < BITS; i++)
//             directionNumbers[1][i] = directionNumbers[1][i-1] ^ (directionNumbers[1][i-1] >> 1);
//     }

// public:
//     SobolSequence(int dim_) : dim(dim_), index(0), state(dim_, 0) {
//         if (dim > 2)
//             throw std::runtime_error(
//                 "Basic Sobol only supports dim<=2. "
//                 "Load Joe-Kuo direction numbers for higher dimensions.");
//         initDirectionNumbers();
//     }

//     // Returns next point in [0,1]^dim
//     std::vector<double> next() {
//         std::vector<double> point(dim);

//         // Gray code increment
//         // int c = __builtin_ctz(index + 1); 

//         int c = 0;
//         int value = index + 1;
//         while ((value & 1) == 0) {
//             value >>= 1;
//             c++;
//         }

//         for (int d = 0; d < dim; d++) {
//             state[d] ^= directionNumbers[d][c];
//             point[d]  = state[d] / std::pow(2.0, BITS);
//         }
//         index++;
//         return point;
//     }

//     void reset() { index = 0; std::fill(state.begin(), state.end(), 0); }
// };



// src/core/SobolSequence.hpp
// #pragma once
// #include <vector>
// #include <cmath>
// #include <stdexcept>

// class SobolSequence {
// private:
//     static const int BITS = 32;
//     int dim;
//     int index;
//     std::vector<unsigned int> state;
//     std::vector<std::vector<unsigned int>> V;  // direction numbers [dim][BITS]

//     void init() {
//         V.resize(dim, std::vector<unsigned int>(BITS, 0));

//         // Dimension 0: standard Van der Corput
//         for (int i = 0; i < BITS; i++)
//             V[0][i] = 1u << (BITS - 1 - i);

//         // Dimensions 1..dim-1: Joe-Kuo primitive polynomials
//         // Each row: s, a, m[1..s]
//         // This is a subset of the full table — extend as needed
//         static const int data[][10] = {
//             // s  a  m1 m2 m3 m4 m5 m6 m7
//             {1, 0, 1},
//             {2, 1, 1, 1},
//             {3, 1, 1, 3, 7},
//             {3, 2, 1, 1, 5},
//             {4, 1, 1, 1, 1, 3},
//             {4, 4, 1, 3, 5, 13},
//             {5, 2, 1, 1, 5, 5, 17},
//             {5, 4, 1, 1, 5, 5, 5},
//             {5, 7, 1, 1, 7, 11, 19},
//             {5, 11, 1, 1, 5, 1, 1},
//             {5, 13, 1, 1, 1, 3, 11},
//             {5, 14, 1, 3, 5, 5, 31},
//             {6, 1, 1, 1, 3, 3, 1, 7},
//             {6, 13, 1, 1, 1, 9, 23, 37},
//             {6, 16, 1, 3, 3, 5, 3, 21},
//             {6, 19, 1, 1, 3, 15, 21, 21},
//             {6, 22, 1, 3, 1, 13, 23, 13},
//             {6, 25, 1, 1, 7, 11, 27, 3},
//             {7, 1, 1, 1, 1, 3, 7, 13, 43},
//             {7, 4, 1, 3, 3, 9, 7, 35, 89},
//             {7, 7, 1, 1, 3, 13, 5, 19, 69},
//             {7, 11, 1, 1, 7, 13, 25, 11, 69},
//             {7, 13, 1, 3, 5, 11, 7, 11, 41},
//             {7, 14, 1, 3, 1, 3, 13, 39, 105},
//             {7, 19, 1, 3, 5, 5, 27, 27, 7},
//             {7, 22, 1, 1, 3, 15, 29, 37, 119},
//             {7, 25, 1, 3, 1, 15, 29, 15, 5},
//             {8, 4, 1, 3, 7, 7, 21, 43, 1, 119},
//             {8, 13, 1, 1, 3, 11, 27, 35, 103, 363},
//             {8, 19, 1, 3, 5, 5, 21, 39, antipodean, 377},
//         };

//         for (int d = 1; d < dim && d <= 30; d++) {
//             int s = data[d-1][0];
//             int a = data[d-1][1];

//             std::vector<int> m(s);
//             for (int i = 0; i < s; i++)
//                 m[i] = data[d-1][2+i];

//             // Scale initial direction numbers
//             for (int i = 0; i < s && i < BITS; i++)
//                 V[d][i] = (unsigned int)m[i] << (BITS - 1 - i);

//             // Recurrence for remaining direction numbers
//             for (int i = s; i < BITS; i++) {
//                 V[d][i] = V[d][i-s] ^ (V[d][i-s] >> s);
//                 for (int k = 1; k < s; k++)
//                     if ((a >> (s-1-k)) & 1)
//                         V[d][i] ^= V[d][i-k];
//             }
//         }

//         state.assign(dim, 0);
//         index = 0;
//     }

// public:
//     SobolSequence(int dim_) : dim(dim_) {
//         if (dim_ > 31)
//             throw std::runtime_error("Extend direction numbers table for dim > 31");
//         init();
//     }

//     std::vector<double> next() {
//         // int c = (index == 0) ? BITS - 1 : __builtin_ctz(index);
//         int c = 0;
//         int value = index;
//         while((value & 1) == 0)
//         {
//             value >>= 1;
//             c++;
//         }
//         std::vector<double> point(dim);
//         for (int d = 0; d < dim; d++) {
//             state[d] ^= V[d][c];
//             point[d]  = state[d] / std::pow(2.0, BITS);
//         }
//         index++;
//         return point;
//     }

//     void reset() {
//         state.assign(dim, 0);
//         index = 0;
//     }
// };

#pragma once
#include <vector>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include "SobolDirectionLoader.hpp"

class SobolSequence {
private:
    static const int BITS = 32;
    int dim;
    int index;
    std::vector<unsigned int> state;
    std::vector<std::vector<unsigned int>> V;

    // Owen-style randomized scrambling (base-2 nested bit permutations).
    // We implement the special case where, at each node, the only permutation
    // is either identity or swap {0,1}. This is sufficient for base-2.
    bool owenEnabled = false;
    std::uint64_t owenSeed = 0;
    int owenDepth = 12; // number of leading bits to scramble (truncated Owen)

    static std::uint64_t splitmix64(std::uint64_t x) {
        x += 0x9E3779B97F4A7C15ULL;
        x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
        x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
        return x ^ (x >> 31);
    }

    unsigned int owenScrambleState(unsigned int raw, int d) const {
        if (!owenEnabled) return raw;
        if (owenDepth <= 0) return raw;
        if (owenDepth > BITS) return raw;

        // Scramble only the leading owenDepth bits: positions [BITS-owenDepth .. BITS-1]
        // p=0 corresponds to MSB (bit 31), p=owenDepth-1 is the last scrambled leading bit.
        unsigned int out = raw;
        for (int p = 0; p < owenDepth; ++p) {
            int bitPos = (BITS - 1) - p;
            unsigned int inputBit = (raw >> bitPos) & 1U;

            // prefix = first p bits (from MSB) of the *original* raw digit sequence
            // For p=0, prefix is empty => 0.
            unsigned int prefix = (p == 0) ? 0U : (raw >> (BITS - p));

            std::uint64_t h = splitmix64(owenSeed ^ (std::uint64_t)d * 0x9E3779B97F4A7C15ULL
                                         ^ (std::uint64_t)p * 0xBF58476D1CE4E5B9ULL
                                         ^ (std::uint64_t)prefix * 0x94D049BB133111EBULL);
            unsigned int flip = static_cast<unsigned int>(h & 1ULL); // 0 identity, 1 swap

            unsigned int outputBit = flip ? (inputBit ^ 1U) : inputBit;
            if (outputBit != inputBit) {
                out ^= (1U << bitPos);
            }
        }
        return out;
    }

    void initFromTable(const std::vector<SobolDirectionData>& table) {
        V.resize(dim, std::vector<unsigned int>(BITS, 0));

        // Dimension 0: Van der Corput
        for (int i = 0; i < BITS; i++)
            V[0][i] = 1u << (BITS - 1 - i);

        // Dimensions 1..dim-1 from loaded table
        for (int d = 1; d < dim; d++) {
            if (d - 1 >= (int)table.size())
                throw std::runtime_error("Not enough dimensions in direction numbers file");

            const auto& row = table[d - 1];
            int s = row.s;
            int a = row.a;

            // Scale initial direction numbers
            for (int i = 0; i < s && i < BITS; i++)
                V[d][i] = (unsigned int)row.m[i] << (BITS - 1 - i);

            // Recurrence for remaining
            for (int i = s; i < BITS; i++) {
                V[d][i] = V[d][i-s] ^ (V[d][i-s] >> s);
                for (int k = 1; k < s; k++)
                    if ((a >> (s - 1 - k)) & 1)
                        V[d][i] ^= V[d][i-k];
            }
        }

        state.assign(dim, 0);
        index = 0;
    }

public:
    // Load from Joe-Kuo file
    SobolSequence(int dim_, const std::string& directionFile)
        : dim(dim_)
    {
        auto table = SobolDirectionLoader::load(directionFile);
        if ((int)table.size() < dim - 1)
            throw std::runtime_error("Direction numbers file has fewer dimensions than requested");
        initFromTable(table);
    }

    // Same as above, but adds Owen-style scrambling (truncated to owenDepth leading bits).
    SobolSequence(int dim_, const std::string& directionFile,
                   std::uint64_t owenSeed_, int owenDepth_)
        : dim(dim_), owenEnabled(true), owenSeed(owenSeed_), owenDepth(owenDepth_)
    {
        auto table = SobolDirectionLoader::load(directionFile);
        if ((int)table.size() < dim - 1)
            throw std::runtime_error("Direction numbers file has fewer dimensions than requested");
        initFromTable(table);
    }

    std::vector<double> next() {
        // Gray code: find rightmost zero bit of index
        // Special case: index=0 means c=0 (first point)
        // int c = 0;
        // if (index > 0) {
        //     int value = index;
        //     while ((value & 1) == 1) {  // find rightmost ZERO bit, not one
        //         value >>= 1;
        //         c++;
        //     }
        // }

        int value = index + 1;
        int c = 0;

        // while ((value & 1) == 0) {
        //     value >>= 1;
        //     c++;
        // }
        
        if (index > 0) {
            int value = index;
            while (value & 1) {
                value >>= 1;
                c++;
            }
        }

        std::vector<double> point(dim);
        for (int d = 0; d < dim; d++) {
            state[d] ^= V[d][c];
            unsigned int raw = state[d];
            unsigned int scr = owenEnabled ? owenScrambleState(raw, d) : raw;
            // point[d]  = state[d] / std::pow(2.0, BITS);
            static const double SCALE = 1.0 / (1ULL << BITS);
            point[d] = (scr + 0.5) * SCALE;
            // point[d] = (state[d] + 0.5) / std::pow(2.0, BITS);
        }
        index++;
        return point;
    }

    void reset() {
        state.assign(dim, 0);
        index = 0;
    }
};