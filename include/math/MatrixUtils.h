#pragma once

#include <tuple>
#include <utility/Concepts.h>

template<numeric T>
class Matrix;

template<numeric T>
class Vector3;

namespace matrix {
    /**
     * @brief Generate a 3x3 extrinsic rotation matrix.
     */
    Matrix<double> rotation_matrix(double alpha, double beta, double gamma);

    /**
     * @brief Generate a 3x3 rotation matrix from a rotation axis and an angle around this axis. 
     *        This uses the Euler-Rodrigues formulation.
     * @param axis The rotation axis.
     * @param angle The rotation angle.
     */
    Matrix<double> rotation_matrix(const Vector3<double>& axis, double angle);

    /**
     * @brief Get the identity matrix of a given dimension. 
     */
    Matrix<double> identity(unsigned int dim);
}

namespace vector3 {
    /**
     * @brief Generate a complete 3D basis from a single basis vector. 
     *        Implementation based on the "frisvad" algorithm from https://backend.orbit.dtu.dk/ws/portalfiles/portal/126824972/onb_frisvad_jgt2012_v2.pdf.
     * 
     * @param v The first basis vector. 
     */
    std::tuple<Vector3<double>, Vector3<double>, Vector3<double>> generate_basis(const Vector3<double>& v);
}