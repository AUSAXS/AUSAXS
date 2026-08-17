// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/Molecule.h>
#include <rigidbody/RigidbodyFwd.h>
#include <io/IOFwd.h>

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace ausaxs::rigidbody::continuation {
    /**
     * @brief The file extension identifying a continuation state file.
     *
     * Written by a `save <name>.continue` element and consumed by `load {continue <name>.continue}`.
     */
    inline constexpr std::string_view continuation_extension = ".continue";

    /**
     * @brief One distance constraint, described independently of the objects it was built from.
     *
     * Constraints hold raw body/atom indices into the molecule they were created against, so they cannot be copied
     * across as objects. Recording what they *mean* instead lets the reader rebuild them against the restored molecule.
     */
    struct ContinuationConstraint {
        enum class Kind : std::uint8_t {bond, cm, attractor, repeller, atom};

        Kind kind = Kind::bond;
        int ibody1 = -1, iatom1 = -1;
        int ibody2 = -1, iatom2 = -1;
        std::pair<int, int> isym1 = {-1, -1};
        std::pair<int, int> isym2 = {-1, -1};
        double d_target = 0;
    };

    /**
     * @brief Everything needed to resume a refinement exactly where a previous one stopped.
     *
     * This is deliberately richer than a .pdb: writing a molecule to .pdb flattens every symmetry into explicit atoms
     * and discards the residue numbering and Cα flags, none of which can be recovered on re-reading. A continuation
     * state instead keeps the bodies as they were — implicit symmetries, per-atom metadata, body names and all — so a
     * follow-up run refines the same system rather than a lookalike.
     */
    struct ContinuationState {
        ContinuationState() = default;
        ContinuationState(ContinuationState&&) = default;
        ContinuationState& operator=(ContinuationState&&) = default;

        data::Molecule molecule;
        std::vector<std::string> body_names;            //< display name per body, parallel to the molecule's bodies
        std::vector<ContinuationConstraint> constraints;
    };

    /**
     * @brief Write the current state of a rigid body to a continuation file.
     *
     * @param path The destination. Its directory is created if missing.
     * @param rigidbody The rigid body whose molecule and constraints are recorded.
     * @param body_names Display name per body, ordered by body index. May be empty, in which case names are not stored.
     */
    void write_continuation_state(const ausaxs::io::File& path, const Rigidbody& rigidbody, const std::vector<std::string>& body_names);

    /**
     * @brief Read a continuation file previously written by write_continuation_state.
     *
     * @throws std::runtime_error if the file is missing, is not a continuation file, or was written by an incompatible version.
     */
    [[nodiscard]] ContinuationState read_continuation_state(const ausaxs::io::File& path);
}
