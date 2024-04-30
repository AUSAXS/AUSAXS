/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/parameters/RotationsOnly.h>
#include <math/Vector3.h>

#include <random>

using namespace rigidbody::parameter;

RotationsOnly::~RotationsOnly() = default;

Parameter RotationsOnly::next() {
    double scaling = decay_strategy->next();
    double dr1 = rotation_dist(generator)*scaling;
    double dr2 = rotation_dist(generator)*scaling;
    double dr3 = rotation_dist(generator)*scaling;
    return Parameter({0, 0, 0}, dr1, dr2, dr3);
}