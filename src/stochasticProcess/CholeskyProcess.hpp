#pragma once
#include "MultiDimensionalProcess.hpp"

class CholeskyProcess : public MultiDimensionalProcess {
public:
    virtual bool diffusionIsCholesky() const final { return true; }
};