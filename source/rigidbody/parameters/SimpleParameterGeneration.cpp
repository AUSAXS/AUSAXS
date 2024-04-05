/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/parameters/SimpleParameterGeneration.h>
#include <math/Vector3.h>

#include <random>

using namespace rigidbody::parameter;

SimpleParameterGeneration::~SimpleParameterGeneration() = default;

Parameter SimpleParameterGeneration::next() {
    double scaling = decay_strategy->next();

    double dx = translation_dist(generator)*scaling;
    double dy = translation_dist(generator)*scaling;
    double dz = translation_dist(generator)*scaling;

    double dr1 = rotation_dist(generator)*scaling;
    double dr2 = rotation_dist(generator)*scaling;
    double dr3 = rotation_dist(generator)*scaling;

    return Parameter(Vector3(dx, dy, dz), dr1, dr2, dr3);
}