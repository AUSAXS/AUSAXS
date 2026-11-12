// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/scattering_model/DistanceHistogram.h>
#include <utility/TypeTraits.h>

namespace ausaxs::hist {
    /**
     * @brief A scattering model where all scatterers are point-like, i.e. they have no form factor and contribute only their weight to the scattering.
     */
    class ScatteringModelPoint : public DistanceHistogram {
        public:
            using DistanceHistogram::DistanceHistogram;

            /**
             * @brief Reinterpret the total histogram of any other histogram or scattering model as point scatterers.
             *        Its distance axis is kept, so weighted bins remain weighted. Only the DistanceHistogram part of @a hist is consumed.
             */
            explicit ScatteringModelPoint(DistanceHistogram&& hist) noexcept;
    };
    static_assert(supports_nothrow_move_v<ScatteringModelPoint>, "ScatteringModelPoint should support nothrow move semantics.");
}
