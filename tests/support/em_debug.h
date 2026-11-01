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
     *
     * That constructor is protected because a stack built without a map file has no header, and every production entry point
     * (cli_em_fitter, em_fitter_gui) reads a file. Tests that want a synthetic stack use this subclass instead, which keeps the
     * unheadered path visible to the tests without offering it as public API.
     *
     * Note that ImageStackBase::get_header() throws on such a stack until set_header() has been called.
     */
    struct ImageStackBaseDebug : public em::ImageStackBase {
        // an inherited constructor keeps the base's access, so it has to be forwarded explicitly
        ImageStackBaseDebug(const std::vector<em::Image>& images) : em::ImageStackBase(images) {}
    };

    /**
     * @brief Test-only access to ImageStack's list-of-images constructor. See ImageStackBaseDebug.
     */
    struct ImageStackDebug : public em::ImageStack {
        ImageStackDebug(const std::vector<em::Image>& images) : em::ImageStack(images) {}
    };
}
