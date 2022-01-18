#pragma once

#include <memory>

#include "data/Protein.h"
#include "rigidbody/Constraint.h"
#include "rigidbody/BodySelectStrategy.h"

class RigidBody {
  public:
    /**
     * @brief Construtor. 
     * 
     * Prepare a new rigid body for optimization. 
     * 
     * @param protein The protein to be optimized. 
     */
    RigidBody(Protein& protein) : protein(protein) {}

    /**
     * @brief Perform a rigid-body optimization for this structure. 
     */
    void optimize();

    /**
     * @brief Add a constraint to this rigid body. 
     */
    void add_constraint(const Constraint& constraint);

    /**
     * @brief Create a constraint for this rigid body. 
     */
    void create_constraint(const std::shared_ptr<Atom> const atom1, const std::shared_ptr<Atom> const atom2);

    Protein& protein;
    std::vector<Constraint> constraints;

  private:
    std::unique_ptr<BodySelectStrategy> body_selector;

    /**
     * @brief Perform a single step of the optimization, and calculate the resulting chi2 value. 
     */
    double chi2();

    /**
     * @brief Rotate a body with the currently chosen transformation strategy. 
     */
    void rotate();

    /**
     * @brief Translate a body with the currently chosen transformation strategy. 
     */
    void translate();
};