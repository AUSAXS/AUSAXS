// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace ausaxs::rigidbody::sequencer {
    // elements accepting an open-ended list of inline arguments declare this as their maximum.
    constexpr std::size_t unbounded_inline_args = 100;

    /**
     * @brief The inline (positional) arguments an element accepts.
     *
     * The names are what the element calls its arguments, in script order. They carry no meaning for the parser - the
     * element still reads its arguments by index - but they are what the error messages and the GUI show the user, so
     * a slot absorbing the remaining arguments should say so, e.g. "others...".
     */
    struct InlineSignature {
        std::vector<std::string> names;
        std::size_t min = 0;
        std::size_t max = 0;
    };
}
