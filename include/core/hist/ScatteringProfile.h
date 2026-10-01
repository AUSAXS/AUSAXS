// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <dataset/DatasetFwd.h>
#include <math/Vector.h>
#include <utility/Axis.h>
#include <utility/TypeTraits.h>

#include <string>
#include <vector>

namespace ausaxs::hist {
    /**
     * @brief A calculated scattering intensity I(q), sampled on a q-axis.
     */
    class ScatteringProfile {
        public:
            ScatteringProfile() = default;
            ScatteringProfile(std::vector<double> I, const Axis& q_axis);

            /**
             * @brief Get the q-axis the intensity is sampled on.
             */
            [[nodiscard]] const Axis& get_axis() const;

            /**
             * @brief Get the intensity values.
             */
            [[nodiscard]] const std::vector<double>& get_intensity() const;

            [[nodiscard]] int size() const noexcept;

            /**
             * @brief Get this profile as a [q | I] dataset.
             */
            [[nodiscard]] Dataset as_dataset() const;

            [[nodiscard]] std::string to_string() const;

            ScatteringProfile& operator+=(const ScatteringProfile& rhs);
            ScatteringProfile& operator-=(const ScatteringProfile& rhs);
            ScatteringProfile& operator*=(double rhs);
            double& operator[](int i);
            double operator[](int i) const;
            bool operator==(const ScatteringProfile& rhs) const;

        private:
            Vector<double> I;
            Axis axis;
    };

    ScatteringProfile operator+(const ScatteringProfile& lhs, const ScatteringProfile& rhs);
    ScatteringProfile operator-(const ScatteringProfile& lhs, const ScatteringProfile& rhs);
    ScatteringProfile operator*(const ScatteringProfile& lhs, double rhs);

    static_assert(supports_nothrow_move_v<ScatteringProfile>, "ScatteringProfile should be noexcept move constructible.");
}
