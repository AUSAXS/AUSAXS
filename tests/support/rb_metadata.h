// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/Molecule.h>
#include <data/Body.h>
#include <data/atoms/AtomMetadata.h>
#include <form_factor/FormFactorType.h>
#include <constants/Constants.h>

namespace ausaxs::test {
    /**
     * @brief Attach C-alpha backbone metadata to every body of an in-memory molecule.
     *
     * The rigidbody constraint code identifies the backbone atoms it constrains through the
     * per-atom metadata that PDB loading attaches when settings::molecule::store_calpha is enabled.
     * Molecules assembled directly from AtomFF vectors in tests never pass through that path, so the
     * metadata is absent and the constraint generators dereference an empty optional.
     *
     * This marks every carbon atom as a C-alpha backbone atom, reproducing the carbon-based atom
     * selection the constraints used before backbone metadata was introduced.
     */
    inline void mark_backbone_carbons(data::Molecule& mol) {
        for (unsigned int b = 0; b < mol.size_body(); ++b) {
            data::Body& body = mol.get_body(b);
            data::AtomMetadata md;
            auto& backbone = md.backbone.emplace();
            backbone.reserve(body.size_atom());
            for (unsigned int i = 0; i < body.size_atom(); ++i) {
                backbone.push_back(
                    form_factor::to_atom_type(body.get_atom(i).form_factor_type()) == constants::atom_t::C
                        ? data::backbone_t::c_alpha
                        : data::backbone_t::none
                );
            }
            body.set_metadata(std::move(md));
        }
    }
}
