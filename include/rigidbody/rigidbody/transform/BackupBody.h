// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/Body.h>
#include <rigidbody/parameters/BodyTransformParametersAbsolute.h>

namespace ausaxs::rigidbody::transform {
    struct BackupBody {
        BackupBody(const data::Body& body, unsigned int index, const parameter::BodyTransformParametersAbsolute& params) 
            : index(index), body(body), params(params) 
        {}

        unsigned int index;
        data::Body body;
        parameter::BodyTransformParametersAbsolute params;
    };
}