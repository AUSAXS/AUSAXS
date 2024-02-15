#pragma once

#include <math/Vector3.h>

namespace rigidbody::parameter {
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
         * @param dr The translation vector.
         * @param alpha The first Euler angle.
         * @param beta The second Euler angle.
         * @param gamma The third Euler angle.
         */
        Parameter(const Vector3<double>& dr, double alpha, double beta, double gamma);

        /**
         * @brief Get a string representation of this Parameter.
         */
        std::string to_string() const;

        /**
         * @brief Output the string representation of this Parameter to a stream.
         */
        friend std::ostream& operator<<(std::ostream& os, const Parameter& p) {os << p.to_string(); return os;}

        Vector3<double> dr;
        double alpha, beta, gamma;
    };
}