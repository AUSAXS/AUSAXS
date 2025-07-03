// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <utility/UtilityFwd.h>

#include <vector>

namespace ausaxs::em {
    /**
     * @brief Describes the bounds of some object contained within a 2D matrix. 
     */
    class ObjectBounds2D {
        public:
            ObjectBounds2D(unsigned int size_x, unsigned int size_y);

            ~ObjectBounds2D();

            /**
             * @brief Set the minimum bound of the xth row.
             */
            void set_min(unsigned int x, unsigned int min);

            /**
             * @brief Set the maximum bound of the xth row.
             */
            void set_max(unsigned int x, unsigned int max);

            /**
             * @brief Set the bounds of the xth row. 
             */
            void set_bounds(unsigned int x, const Limit& limit);

            /**
             * @brief Set the bounds of the xth row. 
             */
            void set_bounds(unsigned int x, unsigned int min, unsigned int max);

            /**
             * @brief Get the bounds of the xth row. 
             */
            [[nodiscard]] const Limit& operator[](unsigned int x) const;

            /**
             * @brief Get the size in the x-dimension. 
             */
            [[nodiscard]] unsigned int size_x() const;

            /**
             * @brief Get the size in the y-dimension. 
             */
            [[nodiscard]] unsigned int size_y() const;

            /**
             * @brief Returns true if no bounds have been set. 
             */
            [[nodiscard]] bool empty() const;

            /**
             * @brief Get the total bounded area.
             */
            [[nodiscard]] unsigned int bounded_area() const;

            /**
             * @brief Get the total area.
             */
            [[nodiscard]] unsigned int total_area() const;

            [[nodiscard]] bool operator==(const ObjectBounds2D& other) const;

        private:
            std::vector<Limit> bounds;
            unsigned int N, M;
    };
}