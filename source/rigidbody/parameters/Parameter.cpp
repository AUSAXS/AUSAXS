/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/parameters/Parameter.h>

using namespace rigidbody::parameter;

Parameter::Parameter() : dr(0, 0, 0), alpha(0), beta(0), gamma(0) {}

Parameter::Parameter(const Vector3<double>& dr, double alpha, double beta, double gamma) : dr(dr), alpha(alpha), beta(beta), gamma(gamma) {}