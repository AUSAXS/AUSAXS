// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <string_view>

namespace ausaxs::gpu::shader::source {
    // The wgsl sources are embedded at build time by embed_shaders.cmake, so no shader files have to
    // be shipped alongside the binaries. The definitions live in the generated EmbeddedShaders.cpp.
    extern const std::string_view simple_weighted;
    extern const std::string_view simple_unweighted;
}
