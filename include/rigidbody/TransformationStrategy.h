#pragma once

#include "data/Protein.h"
#include "rigidbody/Constraint.h"

// forwards declaration
class RigidBody;

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
    TransformationStrategy(const RigidBody* protein) : protein(protein) {}

    /**
     * @brief Destructor.
     */
    virtual ~TransformationStrategy() = default;

    /**
     * @brief Rotate a body. 
     * 
     * @param angle The rotation angle in radians. 
     * @param body The body being rotated. 
     */
    virtual void rotate(const double rad, Constraint& constraint) = 0;

    /**
     * @brief Translate a body. 
     * 
     * @param length The distance to translate. 
     * @param body The body being translated. 
     */
    virtual void translate(const double& length, Constraint& constraint) = 0;

    /**
     * @brief Get all bodies connected by constraints to the first body of the pivot. 
     *        If we have the four bodies A - B - C - D and pivot around the BC connection, this would return the group {AB}.
     */
    vector<Body*> get_connected(const Constraint& pivot) const;

  protected: 
    const RigidBody* protein;
};