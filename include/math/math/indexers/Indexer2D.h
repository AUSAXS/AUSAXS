// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <cassert>

#include <math/detail/Diagnostics.h>
#include <math/indexers/Shape.h>

namespace ausaxs::utility::indexer {
    /**
     * @brief CRTP mixin providing element access for a two-dimensional container.
     *        The deriving class must expose a contiguous @c data member (row-major, with row length @c M) and the dimensions @c N and @c M. 
     *        With Shape::Triangular, @c N and @c M must be equal and @c data holds one element per unordered pair (i, j).
     */
    template<typename Derived, Shape S = Shape::Square>
    class Indexer2D {
        Indexer2D() = default;
        friend Derived;

        protected:
            constexpr const auto& index(int i, int j) const {
                assert([&]() -> bool {
                    if (0 <= i && i < derived().N && 0 <= j && j < derived().M) {return true;}
                    return ausaxs::detail::report_index_2d("Indexer2D", i, j, derived().N, derived().M);
                }() && "Indexer2D: Index out of bounds.");
                return derived().data[pair_offset(i, j)];
            }

            constexpr auto& index(int i, int j) {
                assert([&]() -> bool {
                    if (0 <= i && i < derived().N && 0 <= j && j < derived().M) {return true;}
                    return ausaxs::detail::report_index_2d("Indexer2D", i, j, derived().N, derived().M);
                }() && "Indexer2D: Index out of bounds.");
                return derived().data[pair_offset(i, j)];
            }

            constexpr const auto& linear_index(int i) const { 
                assert([&]() -> bool {
                    if (0 <= i && i < pair_count()) {return true;}
                    return ausaxs::detail::report_index_1d("Indexer2D::linear_index", i, pair_count());
                }() && "Indexer2D::linear_index: Index out of bounds.");
                return derived().data[i];
            }

            constexpr auto& linear_index(int i) { 
                assert([&]() -> bool {
                    if (0 <= i && i < pair_count()) {return true;}
                    return ausaxs::detail::report_index_1d("Indexer2D::linear_index", i, pair_count());
                }() && "Indexer2D::linear_index: Index out of bounds.");
                return derived().data[i];
            }

            /**
             * @brief The position of the element of the pair (i, j) in @c data.
             */
            constexpr int pair_offset(int i, int j) const {
                if constexpr (S == Shape::Square) {return square::pair_index(i, j, derived().M);}
                assert(derived().N == derived().M && "Indexer2D: a triangular container must have equal dimensions.");
                return triangular::pair_index(i, j, derived().N);
            }

            /**
             * @brief The number of stored elements.
             */
            constexpr int pair_count() const {
                if constexpr (S == Shape::Square) {return square::pair_count(derived().N, derived().M);}
                assert(derived().N == derived().M && "Indexer2D: a triangular container must have equal dimensions.");
                return triangular::pair_count(derived().N);
            }

        private:
            Derived& derived() { return static_cast<Derived&>(*this); }
            const Derived& derived() const { return static_cast<const Derived&>(*this); }
    };
}