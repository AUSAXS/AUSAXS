#pragma once

#include <math/Vector3.h>

#include <unordered_map>
#include <vector>

class Protein;
namespace rigidbody {
    /**
     * @brief A small structure for storing a single set of parameters. 
     */
    struct Parameter {
        /**
         * @brief Default constructor.
         */
        Parameter();

        /**
         * @brief Constructor.
         * 
         * @param dx The translation vector.
         * @param alpha The first Euler angle.
         * @param beta The second Euler angle.
         * @param gamma The third Euler angle.
         */
        Parameter(const Vector3<double>& dx, double alpha, double beta, double gamma);

        /**
         * @brief Get a string representation of this Parameter.
         */
        std::string to_string() const;

        /**
         * @brief Output the string representation of this Parameter to a stream.
         */
        friend std::ostream& operator<<(std::ostream& os, const Parameter& p) {os << p.to_string(); return os;}

        Vector3<double> dx;
        double alpha, beta, gamma;
    };

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