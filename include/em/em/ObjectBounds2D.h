// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <utility/Limit.h>

#include <vector>

namespace ausaxs::em {
    /**
     * @brief Describes the bounds of some object contained within a 2D matrix. 
     *
     * Each row is bounded by a half-open index range [min, max), so @a max is one past the last enclosed index and
     * an empty row is expressed as an empty range. A newly constructed instance encloses the entire matrix.
     */
    class ObjectBounds2D {
        public:
            ObjectBounds2D(int size_x, int size_y);

            ~ObjectBounds2D();

            /**
             * @brief Set the minimum bound of the xth row.
             */
            void set_min(int x, int min);

            /**
             * @brief Set the maximum bound of the xth row, exclusive.
             */
            void set_max(int x, int max);

            /**
             * @brief Set the bounds of the xth row. 
             */
            void set_bounds(int x, const Limit& limit);

            /**
             * @brief Set the bounds of the xth row. 
             */
            void set_bounds(int x, int min, int max);

            /**
             * @brief Get the bounds of the xth row. 
             */
            [[nodiscard]] const Limit& operator[](int x) const;

            /**
             * @brief Get the size in the x-dimension. 
             */
            [[nodiscard]] int size_x() const;

            /**
             * @brief Get the size in the y-dimension. 
             */
            [[nodiscard]] int size_y() const;

            /**
             * @brief Returns true if no area is enclosed by these bounds.
             */
            [[nodiscard]] bool empty() const;

            /**
             * @brief Get the total bounded area.
             */
            [[nodiscard]] int bounded_area() const;

            /**
             * @brief Get the total area.
             */
            [[nodiscard]] int total_area() const;

            [[nodiscard]] bool operator==(const ObjectBounds2D& other) const;

        private:
            std::vector<Limit> bounds;
            int N, M;
    };
}