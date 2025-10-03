// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/detail/Configuration.h>
#include <rigidbody/Rigidbody.h>
#include <data/Body.h>

using namespace ausaxs;

rigidbody::detail::Configuration::Configuration() : chi2(std::numeric_limits<double>::max()) {}

rigidbody::detail::Configuration::Configuration(observer_ptr<const Rigidbody> rigidbody) noexcept 
    : parameters(rigidbody->molecule.size_body()), chi2(std::numeric_limits<double>::max())
{
    for (int i = 0; i < static_cast<int>(parameters.size()); ++i) {
        parameters[i].symmetry_pars.resize(rigidbody->molecule.get_body(i).size_symmetry());
    }
}

rigidbody::detail::Configuration::Configuration(observer_ptr<const Rigidbody> rigidbody, double chi2) noexcept 
    : parameters(rigidbody->molecule.size_body()), chi2(chi2) 
{}

rigidbody::detail::Configuration::~Configuration() = default;