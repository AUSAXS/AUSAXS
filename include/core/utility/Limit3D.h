// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <utility/Limit.h>

namespace ausaxs {
    /**
     * @brief A representation of a 3-dimensional limited span of values. 
     */
    class Limit3D {
        public:
            Limit3D() noexcept;
            /**
             * @brief Construct a new 3D limit from three 1D limits. 
             * 
             * @param x The limit on the x-axis. 
             * @param y The limit on the y-axis. 
             * @param z The limit on the z-axis. 
             */
            Limit3D(const Limit& x, const Limit& y, const Limit& z) noexcept;

            /**
             * @brief Construct a new 3D limit from a minimum and maximum value along each axis. 
             */
            Limit3D(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) noexcept;

            /**
             * @brief Check if this object is fully initialized. Returns false if any of its Limits are empty.
             */
            [[nodiscard]] bool empty() const noexcept;

            Limit x;
            Limit y;
            Limit z;
    };
}