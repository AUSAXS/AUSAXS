// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/Vector3.h>
#include <constants/ConstantsCoordinates.h>
#include <form_factor/FormFactorType.h>

#include <type_traits>

namespace ausaxs::data {
    namespace detail {
        template<typename T>
        struct AtomForwarder {
            using precision_t = constants::coords_precision_t;

            [[nodiscard]] const Vector3<precision_t>& coordinates() const {return static_cast<const T*>(this)->get_atom().coords;}
            [[nodiscard]] Vector3<precision_t>& coordinates() {return static_cast<T*>(this)->get_atom().coords;}
            [[nodiscard]] const Vector3<precision_t>& position() const {return static_cast<const T*>(this)->get_atom().coords;}
            [[nodiscard]] Vector3<precision_t>& position() {return static_cast<T*>(this)->get_atom().coords;}
            [[nodiscard]] precision_t weight() const {return static_cast<const T*>(this)->get_atom().w;}
            [[nodiscard]] precision_t& weight() {return static_cast<T*>(this)->get_atom().w;}

            [[nodiscard]] precision_t x() const {return static_cast<const T*>(this)->get_atom().coords.x();}
            [[nodiscard]] precision_t y() const {return static_cast<const T*>(this)->get_atom().coords.y();}
            [[nodiscard]] precision_t z() const {return static_cast<const T*>(this)->get_atom().coords.z();}
            [[nodiscard]] precision_t& x() {return static_cast<T*>(this)->get_atom().coords.x();}
            [[nodiscard]] precision_t& y() {return static_cast<T*>(this)->get_atom().coords.y();}
            [[nodiscard]] precision_t& z() {return static_cast<T*>(this)->get_atom().coords.z();}

            [[nodiscard]] bool operator==(const AtomForwarder& rhs) const = default;
        };
    }

    /**
     * @brief The most basic information of an atom that is needed to calculate a distance histogram.
     */
    struct Atom : detail::AtomForwarder<Atom> {
        Atom() = default;
        Atom(Vector3<precision_t> coords, precision_t weight) : coords(std::move(coords)), w(weight) {}
        [[nodiscard]] const Atom& get_atom() const {return *this;}
        [[nodiscard]] Atom& get_atom() {return *this;}
        bool operator==(const Atom& rhs) const = default;

        Vector3<precision_t> coords;
        precision_t w;
    };
    static_assert(sizeof(Atom) == 4*sizeof(constants::coords_precision_t), "Atom size is off");
    static_assert(std::is_trivial_v<Atom>,                                 "Atom is not trivial");
    static_assert(std::is_standard_layout_v<Atom>,                         "Atom is not standard layout");
    static_assert(supports_nothrow_move_v<Atom>,                           "Atom should support nothrow move semantics.");

    template<typename T> concept AtomVector = std::is_same_v<std::remove_cvref_t<T>, std::vector<data::Atom>>;
}