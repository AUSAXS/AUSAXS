// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/DataFwd.h>
#include <io/IOFwd.h>

#include <vector>

namespace ausaxs::rigidbody {
    struct BodySplitter {
        /**
         * @brief Load the structural data from a file and split it into multiple bodies at the designated indices.
         */
        static data::Molecule split(const io::File& input, const std::vector<int>& splits);

        /**
         * @brief Load the structural data from a file and split it into multiple bodies based on the chainID. 
         */
        static data::Molecule split(const io::File& input);
    };
}