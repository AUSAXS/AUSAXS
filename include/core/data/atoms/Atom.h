// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <form_factor/FormFactorType.h>
#include <math/Vector3.h>

#include <type_traits>

namespace ausaxs::data {
    namespace detail {
        template<typename T>
        struct AtomForwarder {
            [[nodiscard]] const Vector3<double>& coordinates() const {return static_cast<const T*>(this)->get_atom().coords;}
            [[nodiscard]] Vector3<double>& coordinates()             {return static_cast<T*>(this)->get_atom().coords;}
            [[nodiscard]] const Vector3<double>& position() const    {return coordinates();}
            [[nodiscard]] Vector3<double>& position()                {return coordinates();}
            [[nodiscard]] double weight() const                      {return static_cast<const T*>(this)->get_atom().w;}
            [[nodiscard]] double& weight()                           {return static_cast<T*>(this)->get_atom().w;}

            [[nodiscard]] double x() const                           {return static_cast<const T*>(this)->get_atom().coords.x();}
            [[nodiscard]] double y() const                           {return static_cast<const T*>(this)->get_atom().coords.y();}
            [[nodiscard]] double z() const                           {return static_cast<const T*>(this)->get_atom().coords.z();}
            [[nodiscard]] double& x()                                {return static_cast<T*>(this)->get_atom().coords.x();}
            [[nodiscard]] double& y()                                {return static_cast<T*>(this)->get_atom().coords.y();}
            [[nodiscard]] double& z()                                {return static_cast<T*>(this)->get_atom().coords.z();}

            [[nodiscard]] bool operator==(const AtomForwarder& rhs) const = default;

            friend T;
            private: AtomForwarder() = default;
        };
    }

    /**
     * @brief The most basic information of an atom that is needed to calculate a distance histogram.
     */
    struct Atom : detail::AtomForwarder<Atom> {
        Atom() = default;
        Atom(Vector3<double> coords, double weight) : coords(coords), w(weight) {}
        [[nodiscard]] const Atom& get_atom() const {return *this;}
        [[nodiscard]] Atom& get_atom() {return *this;}
        bool operator==(const Atom& rhs) const = default;

        Vector3<double> coords;
        double w;
    };
    static_assert(sizeof(Atom) == 4*sizeof(double),                        "Atom size is off");
    static_assert(std::is_trivial_v<Atom>,                                 "Atom is not trivial");
    static_assert(std::is_standard_layout_v<Atom>,                         "Atom is not standard layout");
    static_assert(supports_nothrow_move_v<Atom>,                           "Atom should support nothrow move semantics.");

    template<typename T> concept AtomVector = std::is_same_v<std::remove_cvref_t<T>, std::vector<data::Atom>>;
}