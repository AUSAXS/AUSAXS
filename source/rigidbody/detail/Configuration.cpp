/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include "rigidbody/parameters/Parameter.h"
#include <rigidbody/detail/Configuration.h>

using namespace ausaxs;

rigidbody::detail::Configuration::Configuration() : parameters(), chi2(std::numeric_limits<double>::max()) {}

rigidbody::detail::Configuration::Configuration(rigidbody::parameter::Parameter&& pars, double chi2) noexcept 
    : parameters(std::move(pars)), chi2(chi2) 
{}

rigidbody::detail::Configuration::~Configuration() = default;