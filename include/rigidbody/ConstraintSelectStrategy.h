#pragma once

#include "rigidbody/Constraint.h"

/**
 * @brief \class ConstraintSelectStrategy. 
 * 
 * This super-class defines the interface for the body selection strategies for the rigid-body optimization. 
 * More specifically its implementations will decide in which order the bodies will be transformed by the optimization algorithm.
 */
class ConstraintSelectStrategy {
  public:
    /**
     * @brief Construtor. 
     */
    ConstraintSelectStrategy(const vector<Constraint>& constraints) : constraints(constraints) {}

    /**
     * @brief Destructor.
     */
    virtual ~ConstraintSelectStrategy() = default;

    /**
     * @brief Get the index of the next body to be transformed. 
     */
    virtual size_t next() = 0;

  protected: 
    const vector<Constraint>& constraints;
};