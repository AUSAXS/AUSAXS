// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/Body.h>
#include <rigidbody/parameters/BodyTransformParametersAbsolute.h>

namespace ausaxs::rigidbody::transform {
    struct BackupBody {
        BackupBody(const data::Body& body, unsigned int index, const parameter::BodyTransformParametersAbsolute& params) 
            : body(body), index(index), params(params) {}
        data::Body body;
        unsigned int index;
        parameter::BodyTransformParametersAbsolute params;
    };
}