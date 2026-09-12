// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/Body.h>
#include <rigidbody/parameters/BodyTransformParametersAbsolute.h>

namespace ausaxs::rigidbody::transform {
    struct BackupBody {
        BackupBody(data::Body body, int index, parameter::BodyTransformParametersAbsolute params) 
            : index(index), body(std::move(body)), params(std::move(params))
        {}

        int index;
        data::Body body;
        parameter::BodyTransformParametersAbsolute params;
    };
}