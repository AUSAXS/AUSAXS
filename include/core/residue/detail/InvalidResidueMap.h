// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <residue/detail/ResidueMap.h>

namespace ausaxs::residue::detail {
    class InvalidResidueMap : public ResidueMap {
        public:
            InvalidResidueMap() = default;
            ~InvalidResidueMap() override = default;

            double get(const AtomKey& /*key*/) override {return 0;}

            constants::atomic_group_t get_atomic_group(const std::string& /*atom_name*/, constants::atom_t /*atom_type*/) override {return constants::atomic_group_t::unknown;}
    };
}