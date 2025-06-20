// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/parameters/Parameter.h>

namespace ausaxs::rigidbody::detail {
    struct Configuration {
        Configuration();
        Configuration(parameter::Parameter&& parameters, double chi2) noexcept;
        ~Configuration();

        parameter::Parameter parameters;
        double chi2;
    };
}