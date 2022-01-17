#pragma once

#include "data/Protein.h"

/**
 * @brief \class BodySelectStrategy. 
 * 
 * This super-class defines the interface for the body selection strategies for the rigid-body optimization. 
 * More specifically its implementations will decide in which order the bodies will be transformed by the optimization algorithm.
 */
class BodySelectStrategy {
  public:
    /**
     * @brief Construtor. 
     */
    BodySelectStrategy(const Protein& protein) : protein(protein) {}

    /**
     * @brief Destructor.
     */
    virtual ~BodySelectStrategy() = default;

    /**
     * @brief Get the index of the next body to be transformed. 
     */
    virtual size_t next() = 0;

  protected: 
    const Protein& protein; // A reference to the protein to be optimized. We need this to access its constituent bodies. 
};