// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <cassert>

#include <math/detail/Diagnostics.h>
#include <math/indexers/RowRange.h>
#include <math/indexers/Shape.h>

#include <span>

namespace ausaxs::utility::indexer {
    /**
     * @brief CRTP mixin providing element access for a three-dimensional container.
     *        The deriving class must expose a contiguous @c data member (laid out so that @c L is the fastest-varying dimension) and the
     *        dimensions @c N, @c M, and @c L.
     *        With Shape::Triangular, @c N and @c M must be equal and @c data holds one row of length @c L per unordered pair (i, j).
     */
    template<typename Derived, Shape S = Shape::Square>
    class Indexer3D {
        // only the deriving class may construct the mixin
        Indexer3D() = default;
        friend Derived;

        protected:
            constexpr const auto& index(int i, int j, int k) const {
                assert([&]() -> bool {
                    if (0 <= i && i < derived().N && 0 <= j && j < derived().M && 0 <= k && k < derived().L) {return true;}
                    return ausaxs::detail::report_index_3d("Indexer3D", i, j, k, derived().N, derived().M, derived().L);
                }() && "Indexer3D: Index out of bounds.");
                return derived().data[k + derived().L * pair_offset(i, j)];
            }

            constexpr auto& index(int i, int j, int k) {
                assert([&]() -> bool {
                    if (0 <= i && i < derived().N && 0 <= j && j < derived().M && 0 <= k && k < derived().L) {return true;}
                    return ausaxs::detail::report_index_3d("Indexer3D", i, j, k, derived().N, derived().M, derived().L);
                }() && "Indexer3D: Index out of bounds.");
                return derived().data[k + derived().L * pair_offset(i, j)];
            }

            constexpr const auto& linear_index(int i) const {
                assert([&]() -> bool {
                    if (0 <= i && i < pair_count()*derived().L) {return true;}
                    return ausaxs::detail::report_index_1d("Indexer3D::linear_index", i, pair_count()*derived().L);
                }() && "Indexer3D::linear_index: Index out of bounds.");
                return derived().data[i];
            }

            constexpr const auto& linear_index(int ij, int k) const {
                return linear_index(ij * derived().L + k);
            }

            constexpr auto& linear_index(int i) {
                assert([&]() -> bool {
                    if (0 <= i && i < pair_count()*derived().L) {return true;}
                    return ausaxs::detail::report_index_1d("Indexer3D::linear_index", i, pair_count()*derived().L);
                }() && "Indexer3D::linear_index: Index out of bounds.");
                return derived().data[i];
            }

            constexpr auto& linear_index(int ij, int k) {
                return linear_index(ij * derived().L + k);
            }

            /**
             * @brief The row of the pair (i, j): its @c L elements along the last dimension.
             *        For Shape::Triangular, (i, j) and (j, i) name the same row.
             */
            auto row(int i, int j) const {
                assert(check_pair(i, j) && "Indexer3D::row: Index out of bounds.");
                return std::span(derived().data.data() + derived().L*pair_offset(i, j), static_cast<std::size_t>(derived().L));
            }

            auto row(int i, int j) {
                assert(check_pair(i, j) && "Indexer3D::row: Index out of bounds.");
                return std::span(derived().data.data() + derived().L*pair_offset(i, j), static_cast<std::size_t>(derived().L));
            }

            /**
             * @brief Every stored row once, in storage order. For Shape::Triangular that is each unordered pair once, so this
             *        is the loop to use whenever the pair (i, j) itself is not needed.
             */
            auto rows() const {return RowRange(derived().data.data(), pair_count(), derived().L);}
            auto rows() {return RowRange(derived().data.data(), pair_count(), derived().L);}

            /**
             * @brief The position of the row of the pair (i, j) among the stored rows.
             */
            constexpr int pair_offset(int i, int j) const {
                if constexpr (S == Shape::Square) {return square::pair_index(i, j, derived().N);}
                assert(derived().N == derived().M && "Indexer3D: a triangular container must have equal pair dimensions.");
                return triangular::pair_index(i, j, derived().N);
            }

            /**
             * @brief The number of stored rows.
             */
            constexpr int pair_count() const {
                if constexpr (S == Shape::Square) {return square::pair_count(derived().N, derived().M);}
                assert(derived().N == derived().M && "Indexer3D: a triangular container must have equal pair dimensions.");
                return triangular::pair_count(derived().N);
            }

        private:
            bool check_pair(int i, int j) const {
                if (0 <= i && i < derived().N && 0 <= j && j < derived().M) {return true;}
                return ausaxs::detail::report_index_2d("Indexer3D::row", i, j, derived().N, derived().M);
            }

            Derived& derived() { return static_cast<Derived&>(*this); }
            const Derived& derived() const { return static_cast<const Derived&>(*this); }
    };

    template<typename Derived>
    using TriangularIndexer3D = Indexer3D<Derived, Shape::Triangular>;
}
