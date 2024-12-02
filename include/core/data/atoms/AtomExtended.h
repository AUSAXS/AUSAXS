#pragma once

#include <data/atoms/AtomFF.h>

namespace ausaxs::data {
    class AtomExtended : public detail::AtomBasicForwarder<AtomExtended> {
        using precision_t = constants::coords_precision_t;
        public:
            AtomBasic& get_atom_basic() {return atom.get_atom_basic();}
            const AtomBasic& get_atom_basic() const {return atom.get_atom_basic();}

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