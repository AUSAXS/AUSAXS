#pragma once

#include "Matrix.h"
#include "Vector.h"

class LinearSolver {
public:
    LinearSolver() = delete;
    virtual ~LinearSolver(){}

    enum Algorithm {Givens};

    /**
     * @brief Solve a linear equation of the form Ax = b through a QR decomposition.
     */
    virtual Vector solve(const Vector& b) const = 0;

private: 
    Algorithm algorithm;
};