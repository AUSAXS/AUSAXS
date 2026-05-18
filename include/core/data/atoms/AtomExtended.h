// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/atoms/AtomFF.h>

namespace ausaxs::data {
    /**
     * @brief An AtomFF extended with an occupancy field.
     *
     * Used where the fractional occupancy of an atomic site must be retained, e.g. for partially
     * occupied crystallographic structures. The layout is kept trivial and padding-free so it can
     * be used interchangeably with the other atom types in the scattering calculations.
     */
    class AtomExtended : public detail::AtomForwarder<AtomExtended> {
        using precision_t = constants::coords_precision_t;
        public:
            Atom& get_atom() {return atom.get_atom();}
            const Atom& get_atom() const {return atom.get_atom();}

            precision_t get_occupancy() const {return occupancy;}
            void set_occupancy(precision_t v) {occupancy = v;} 

        private:
            AtomFF atom;
            constants::coords_precision_t occupancy;
    };
    static_assert(sizeof(AtomExtended) == sizeof(AtomFF) + sizeof(constants::coords_precision_t), "AtomExtended should have no padding.");
    static_assert(std::is_trivial_v<AtomExtended>,          "AtomExtended is not trivial");
    static_assert(std::is_standard_layout_v<AtomExtended>,  "AtomExtended is not standard layout");
    static_assert(supports_nothrow_move_v<AtomExtended>,    "AtomExtended should support nothrow move semantics.");
}