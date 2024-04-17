/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/parameters/TranslationsOnly.h>
#include <math/Vector3.h>

#include <random>

using namespace rigidbody::parameter;

TranslationsOnly::~TranslationsOnly() = default;

Parameter TranslationsOnly::next() {
    double scaling = decay_strategy->next();

    double dx = translation_dist(generator)*scaling;
    double dy = translation_dist(generator)*scaling;
    double dz = translation_dist(generator)*scaling;

    return Parameter({dx, dy, dz}, 0, 0, 0);
}