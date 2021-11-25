#pragma once

#include <initializer_list>
#include "math/VectorN.h"

class Vector3 : public VectorN {
public:
    using VectorN::VectorN;
    Vector3() : VectorN(3) {}
    Vector3(VectorN v) : VectorN(v) {}

    double& x = data[0];
    double& y = data[1];
    double& z = data[2];
};