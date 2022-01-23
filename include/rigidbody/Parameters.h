#pragma once

#include "math/Vector3.h"
#include "data/Protein.h"
#include <unordered_map>
#include <vector>

using std::vector;

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
        Parameter(const Vector3& dx, const double& rx, const double& ry, const double& rz) : dx(dx), rx(rx), ry(ry), rz(rz) {}
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
     * @brief Update the parameter set for a single body. 
     * 
     * @param uid The unique identifier of the body. 
     * @param param The new set of parameters. 
     */
    void update(unsigned int uid, const Parameter& param);

    /**
     * @brief Get the parameter set for a single body. 
     * 
     * @param uid The unique identifier of the body. 
     */
    const Parameter get(unsigned int uid);

    std::unordered_map<unsigned int, unsigned int> id_to_index;
    vector<Parameter> params;
};