/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/detail/Configuration.h>
#include <rigidbody/Rigidbody.h>

using namespace ausaxs;

rigidbody::detail::Configuration::Configuration() : chi2(std::numeric_limits<double>::max()) {}

rigidbody::detail::Configuration::Configuration(observer_ptr<Rigidbody> rigidbody, double chi2) noexcept 
    : parameters(rigidbody->molecule.size_body()), chi2(chi2) 
{}

rigidbody::detail::Configuration::~Configuration() = default;