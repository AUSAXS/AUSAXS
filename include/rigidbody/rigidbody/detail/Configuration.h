// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/RigidbodyFwd.h>
#include <rigidbody/parameters/BodyTransformParameters.h>
#include <utility/observer_ptr.h>

namespace ausaxs::rigidbody::detail {
    struct Configuration {
        Configuration();
        Configuration(observer_ptr<const Rigidbody> rigidbody) noexcept;
        Configuration(observer_ptr<const Rigidbody> rigidbody, double chi2) noexcept;
        ~Configuration();

        std::vector<parameter::BodyTransformParameters> parameters;
        double chi2;
    };
}