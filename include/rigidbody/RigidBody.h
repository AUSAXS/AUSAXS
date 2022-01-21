#pragma once

#include <memory>

#include "data/Protein.h"
#include "rigidbody/Constraint.h"
#include "rigidbody/BodySelectStrategy.h"
#include "rigidbody/TransformationStrategy.h"
#include "fitter/IntensityFitter.h"

class RigidBody {
  /**
   * @brief \struct Parameters.
   * 
   * A small structure for storing the current set of parameters. 
   */
  struct Parameters {
    /**
     * @brief \struct Parameter. 
     * 
     * A small structure for storing a single set of parameters. 
     */
    struct Parameter {
      Vector3 dx = {0, 0, 0};
      double rx = 0, ry = 0, rz = 0;
    };

    /**
     * @brief Constructor.
     * 
     * Create a new storage container for the parameters. 
     * 
     * @param protein The protein to create this object for. 
     */
    Parameters(const Protein& protein);

    /**
     * @brief Update the parameter set for a single body. 
     * 
     * @param uid The unique identifier of the body. 
     * @param dx The new offset position vector. 
     * @param drx The new offset rotation about the x-axis. 
     * @param dry The new offset rotation about the y-axis. 
     * @param drz The new offset rotation about the z-axis. 
     */
    void update(unsigned int uid, Vector3 dx, double drx, double dry, double drz);

    /**
     * @brief Get the parameter set for a single body. 
     * 
     * @param uid The unique identifier of the body. 
     */
    const Parameter get(unsigned int uid);

    std::unordered_map<unsigned int, unsigned int> id_to_index;
    vector<Parameter> params;
  };

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
     * 
     * This method is linear in the total number of atoms. For constant efficiency, also provide pointers to the bodies the two atoms are part of. 
     */
    void create_constraint(const Atom* const atom1, const Atom* const atom2);

    /**
     * @brief Create a constraint for this rigid body. 
     */
    void create_constraint(const Atom* const atom1, const Atom* const atom2, const Body* const body1, const Body* const body2);

    /**
     * @brief Create a constraint for this rigid body. 
     * 
     * This method is linear in the total number of atoms. For constant efficiency, also provide the bodies the two atoms are part of. 
     */
    void create_constraint(const Atom& atom1, const Atom& atom2);

    /**
     * @brief Create a constraint for this rigid body. 
     */
    void create_constraint(const Atom& atom1, const Atom& atom2, const Body& body1, const Body& body2);

    Protein& protein;
    std::vector<Constraint> constraints;

  private:
    std::unique_ptr<BodySelectStrategy> body_selector;
    std::unique_ptr<TransformationStrategy> transform;

    void driver(const string& measurement_path);

    /**
     * @brief Perform a single step of the optimization, and calculate the resulting chi2 value. 
     */
    double chi2(IntensityFitter& fitter) const;

    /**
     * @brief Rotate a body with the currently chosen transformation strategy. 
     */
    void rotate();

    /**
     * @brief Translate a body with the currently chosen transformation strategy. 
     */
    void translate();

    /**
     * @brief Find the bodies containing the argument atoms. Helper method for create_constraint. 
     *        This is linear in the total number of atoms. 
     * 
     * @param atom1 The first return value will be a pointer to the host body of this atom. 
     * @param atom2 The second return value will be a pointer to the host body of this atom. 
     */
    std::pair<const Body*, const Body*> find_host_bodies(const Atom* const atom1, const Atom* const atom2) const noexcept(false);
};