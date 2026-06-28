// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <cstdint>
#include <optional>
#include <vector>

namespace ausaxs::data {
    /**
     * @brief Backbone classification of an atom.
     */
    enum class backbone_t : std::uint8_t {
        none,       //< not a backbone atom
        c_alpha,    //< C-alpha backbone carbon
    };

    /**
     * @brief Optional per-atom metadata for a Body.
     *        Each field is parallel-indexed to the Body's atom vector and is only engaged if the corresponding setting (see settings::molecule) is enabled at load time.
     */
    struct AtomMetadata {
        std::optional<std::vector<backbone_t>> backbone;     //< engaged iff settings::molecule::store_calpha
        std::optional<std::vector<int>>        residue_seq;  //< residue sequence id; engaged iff settings::molecule::store_calpha (used by structure previews)
        std::optional<std::vector<float>>      occupancy;    //< engaged iff settings::molecule::store_occupancy
    };
}
