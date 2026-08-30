// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <gpu/WebGPU/shaders/ShaderStorage.h>
#include <gpu/WebGPU/shaders/ShaderSource.h>

using namespace ausaxs::gpu;

const ShaderDefinition& shader::Simple::weighted() {
    static ShaderDefinition definition(source::simple_weighted);
    return definition;
}

const ShaderDefinition& shader::Simple::unweighted() {
    static ShaderDefinition definition(source::simple_unweighted);
    return definition;
}
