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
         * Each split id marks the *first* residue of a new body, and is consumed the first time it is encountered,
         * so a residue id repeated later in the body (e.g. a second chain) does not trigger a second split. Every id
         * must therefore name a residue that actually occurs, and must not be the lowest one; anything that would
         * produce an empty body - or that falls in a gap in the residue numbering - is an error rather than a body.
         *
         * The body's per-atom metadata is sliced in lockstep with the atoms. Any explicit hydration on the body is
         * dropped, as there is no well-defined assignment of waters to the resulting bodies.
         *
         * @param body   The body to split. Must carry residue sequence metadata (see settings::molecule::store_residue_seq).
         * @param splits Residue sequence ids to split at; produces splits.size()+1 bodies.
         */
        static std::vector<data::Body> split(const data::Body& body, const std::vector<int>& splits);

        /**
         * @brief Load the structural data from a file and split it into multiple bodies at the designated residue sequence ids.
         *
         * Equivalent to loading the file as a single body and calling the overload above on it.
         */
        static data::Molecule split(const io::File& input, const std::vector<int>& splits);

        /**
         * @brief Load the structural data from a file and split it into multiple bodies based on the chainID.
         *
         * This has no runtime counterpart, since chain identifiers are not retained past load time - and need not be,
         * as distinct chains are exactly what this produces distinct bodies for.
         */
        static data::Molecule split(const io::File& input);
    };
}
