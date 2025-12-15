// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/atoms/Atom.h>
#include <form_factor/FormFactor.h>

namespace ausaxs::data {
    /**
     * @brief The most basic information of an atom that is needed to calculate a distance histogram.
     */
    struct Water : detail::AtomForwarder<Water> {
        Water() = default;
        Water(Vector3<precision_t> coords) : 
            coords(std::move(coords)), 
            w(form_factor::lookup::atomic::raw::get(static_cast<form_factor::form_factor_t>(form_factor::water_bin)).I0()) 
        {}

        form_factor::form_factor_t form_factor_type() const {return form_factor::form_factor_t::OH;}
        [[nodiscard]] const Water& get_atom() const {return *this;}
        [[nodiscard]] Water& get_atom() {return *this;}
        bool operator==(const Water& rhs) const = default;

        Vector3<precision_t> coords;
        precision_t w;
    };
    static_assert(sizeof(Water) == 4*sizeof(constants::coords_precision_t), "WaterBasic size is off");
    static_assert(std::is_trivial_v<Water>,                                 "WaterBasic is not trivial");
    static_assert(std::is_standard_layout_v<Water>,                         "WaterBasic is not standard layout");
    static_assert(supports_nothrow_move_v<Water>,                           "WaterBasic should support nothrow move semantics.");

	template<typename T>
	concept WaterVector = std::is_same_v<std::remove_cvref_t<T>, std::vector<data::Water>>;
}