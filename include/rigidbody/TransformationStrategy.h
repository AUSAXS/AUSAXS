#pragma once

#include "data/Protein.h"

/**
 * @brief \class TransformationStrategy. 
 * 
 * This super-class defines the interface for the body transformation strategies for the rigid-body optimization. 
 * More specifically its implementations essentially specifies how other connected bodies are affected by a transformation. 
 */
class TransformationStrategy {
  public:
    /**
     * @brief Construtor. 
     */
    TransformationStrategy(const Protein& protein) : protein(protein) {}

    /**
     * @brief Destructor.
     */
    virtual ~TransformationStrategy() = default;

    /**
     * @brief Rotate a body. 
     * 
     * @param axis The rotation axis. 
     * @param angle The rotation angle in radians. 
     * @param body The body being rotated. 
     */
    virtual void rotate(const Vector3& axis, const double rad, Body& body) = 0;

    /**
     * @brief Translate a body. 
     * 
     * @param v The translation vector. 
     * @param body The body being translated. 
     */
    virtual void translate(const Vector3& v, Body& body) = 0;

  protected: 
    const Protein& protein; // A reference to the protein to be optimized. We need this to access its constituent bodies. 
};