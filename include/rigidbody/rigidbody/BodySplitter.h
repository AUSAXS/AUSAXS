// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/DataFwd.h>
#include <io/IOFwd.h>

#include <vector>

namespace ausaxs::rigidbody {
    struct BodySplitter {
        /**
         * @brief Split a body into multiple bodies at the designated residue sequence ids.
         *
         * Each split id marks the *first* residue of a new body, and is consumed the first time it is encountered, so a residue id repeated later in the body 
         * (e.g. a second chain) does not trigger a second split. Any explicit waters are dropped. 
         *
         * @param body   The body to split. Must carry residue sequence metadata.
         * @param splits Residue sequence ids to split at; produces splits.size()+1 bodies.
         */
        static std::vector<data::Body> split(const data::Body& body, const std::vector<int>& splits);

        /**
         * @brief Load the structural data from a file and split it into multiple bodies at the designated residue sequence ids.
         */
        static data::Molecule split(const io::File& input, const std::vector<int>& splits);

        /**
         * @brief Load the structural data from a file and split it into multiple bodies based on the chainID.
         */
        static data::Molecule split(const io::File& input);
    };
}
