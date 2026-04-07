// src/barriers/Barrier.hpp
#pragma once

enum class BarrierType {
    UpAndOut,    // knocked out if S crosses H from below
    UpAndIn,     // activated only if S crosses H from below
    DownAndOut,  // knocked out if S crosses H from above
    DownAndIn    // activated only if S crosses H from above
};

struct Barrier {
    double H;           // barrier level
    BarrierType type;
    double rebate;      // paid if knocked out (often 0)

    Barrier(double H_, BarrierType type_, double rebate_ = 0.0)
        : H(H_), type(type_), rebate(rebate_) {}

    bool isOut() const {
        return type == BarrierType::UpAndOut ||
               type == BarrierType::DownAndOut;
    }

    bool isUp() const {
        return type == BarrierType::UpAndOut ||
               type == BarrierType::UpAndIn;
    }
};