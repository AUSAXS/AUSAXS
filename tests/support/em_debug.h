// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <em/ImageStack.h>
#include <em/detail/ImageStackBase.h>
#include <em/Image.h>

#include <vector>

namespace ausaxs::test {
    /**
     * @brief Test-only access to ImageStackBase's list-of-images constructor.
     */
    struct ImageStackBaseDebug : public em::ImageStackBase {
        ImageStackBaseDebug(const std::vector<em::Image>& images) : em::ImageStackBase(images) {}
    };

    /**
     * @brief Test-only access to ImageStack's list-of-images constructor.
     */
    struct ImageStackDebug : public em::ImageStack {
        ImageStackDebug(const std::vector<em::Image>& images) : em::ImageStack(images) {}
    };
}
