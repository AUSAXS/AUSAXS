#pragma once

#include <data/atoms/Atom.h>
#include <data/atoms/AtomFF.h>
#include <data/atoms/Water.h>
#include <type_traits>

namespace ausaxs::data {
    namespace detail {
        template<typename T>
        using StrippedType = std::remove_cvref_t<std::remove_pointer_t<T>>;
    }

    template<typename T>
    concept AtomType = std::is_base_of_v<
        detail::AtomForwarder<detail::StrippedType<T>>, 
        detail::StrippedType<T>
    >;
}