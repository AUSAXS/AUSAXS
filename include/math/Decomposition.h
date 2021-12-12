#pragma once

#include "Matrix.h"
#include "Vector.h"

class Decomposition {
public: 
    Decomposition() {}
    virtual ~Decomposition() {}
    virtual void decompose() = 0;

protected:
    static constexpr double precision = 1e-9;
};