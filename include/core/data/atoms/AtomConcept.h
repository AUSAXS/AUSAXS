#pragma once

#include <data/atoms/AtomFF.h>

#include <concepts>

namespace ausaxs {
    template<typename C> concept concept_atom_t = requires(C a) {
        std::is_same_v<C, data::AtomFF> || std::is_same_v<C, data::Atom> || std::is_same_v<C, data::Water>;
    };
}