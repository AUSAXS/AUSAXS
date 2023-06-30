#pragma once

#include <math/Vector3.h>
#include <rigidbody/parameters/Parameter.h>

#include <unordered_map>
#include <vector>

class Protein;
namespace rigidbody {
    /**
     * @brief A small structure for storing the current set of parameters. 
     */
    struct Parameters {
        /**
         * @brief Constructor.
         * 
         * Create a new storage container for the parameters. 
         * 
         * @param protein The protein to create this object for. 
         */
        Parameters(const Protein* protein);

        /**
         * @brief Update the parameter set for a single body. 
         * 
         * @param uid The unique identifier of the body. 
         * @param dx The new offset position vector. 
         * @param drx The new offset rotation about the x-axis. 
         * @param dry The new offset rotation about the y-axis. 
         * @param drz The new offset rotation about the z-axis. 
         */
        void update(unsigned int uid, Vector3<double> dx, double drx, double dry, double drz);

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
        std::vector<Parameter> params;
    };
}