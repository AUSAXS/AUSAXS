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
        BackupBody(data::Body&& body, unsigned int index, parameter::BodyTransformParametersAbsolute&& params) 
            : index(index), body(std::move(body)), params(std::move(params))
        {}

        unsigned int index;
        std::optional<data::Body> body;
        parameter::BodyTransformParametersAbsolute params;
    };
}