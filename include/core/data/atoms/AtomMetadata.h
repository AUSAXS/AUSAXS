// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <cstdint>
#include <optional>
#include <vector>

namespace ausaxs::data {
    /**
     * @brief Backbone classification of an atom, used to decide whether it is a valid constraint target.
     */
    enum class backbone_t : std::uint8_t {
        none,       //< not a backbone atom
        c_alpha,    //< C-alpha backbone carbon, valid for constraining
    };

    /**
     * @brief Optional per-atom metadata for a Body.
     *
     * Each field is parallel-indexed to the Body's atom vector and is only engaged if the
     * corresponding setting (see settings::molecule) was enabled at load time. This keeps
     * feature-specific information off the hot-path AtomFF type while still living next to
     * the atoms it describes. The owning optional on Body must be cleared whenever atoms are
     * added or removed, as the parallel indices would otherwise desynchronize.
     */
    struct AtomMetadata {
        std::optional<std::vector<backbone_t>> backbone;    //< engaged iff settings::molecule::store_calpha
        std::optional<std::vector<float>>      occupancy;   //< engaged iff settings::molecule::store_occupancy
    };
}
