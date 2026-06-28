// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/Body.h>
#include <data/atoms/AtomMetaData.h>
#include <form_factor/FormFactorType.h>
#include <constants/Constants.h>

namespace ausaxs::rigidbody::constraints {
    /**
     * @brief Decide whether atom @p i of @p body is a valid target for an inter-body distance constraint.
     *
     * Prefers the C-alpha backbone metadata when it is available (see settings::molecule::store_calpha),
     * which precisely identifies backbone C-alpha atoms. When no metadata is present it falls back to the
     * legacy heuristic of selecting any carbon atom.
     */
    inline bool is_constraint_candidate(const data::Body& body, unsigned int i) {
        const auto& md = body.get_metadata();
        if (md && md->backbone) {
            return (*md->backbone)[i] == data::backbone_t::c_alpha;
        }
        return form_factor::to_atom_type(body.get_atom(i).form_factor_type()) == constants::atom_t::C;
    }
}
