// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/scattering_model/DistanceHistogram.h>
#include <utility/TypeTraits.h>

namespace ausaxs::hist {
    /**
     * @brief Point scatterers: the plain Debye sum of the total histogram, without the exp(-q²) form factor of DistanceHistogram.
     *        Each scatterer contributes only its weight, so this holds at any q.
     *        Like DistanceHistogram, it keeps only the total histogram, so it is not an IScatteringModel.
     */
    class ScatteringModelPoint : public DistanceHistogram {
        public:
            using DistanceHistogram::DistanceHistogram;

            /**
             * @brief Reinterpret the total histogram of any other histogram or scattering model as point scatterers.
             *        Its distance axis is kept, so weighted bins remain weighted. Only the DistanceHistogram part of @a hist is consumed.
             */
            explicit ScatteringModelPoint(DistanceHistogram&& hist) noexcept;

            ScatteringProfile debye_transform() const override;
            using DistanceHistogram::debye_transform;

        protected:
            std::vector<double> debye_transform_beyond_axis(const std::vector<double>& q) const override;
    };
    static_assert(supports_nothrow_move_v<ScatteringModelPoint>, "ScatteringModelPoint should support nothrow move semantics.");
}
