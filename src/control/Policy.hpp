// src/control/Policy.hpp
#pragma once

class Policy {
public:
    virtual double control(double x, double t) const = 0;
    virtual ~Policy() {}
};