// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <cassert>

#include <math/detail/Diagnostics.h>

namespace ausaxs::utility::indexer {
    /**
     * @brief CRTP mixin providing element access for a three-dimensional container.
     *        The deriving class must expose a contiguous @c data member (laid out so that @c L is the fastest-varying dimension) and the 
     *        dimensions @c N, @c M, and @c L. 
     */
    template<typename Derived>
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
                return derived().data[k + derived().L * (j + derived().M * i)]; 
            }

            constexpr auto& index(int i, int j, int k) {
                assert([&]() -> bool {
                    if (0 <= i && i < derived().N && 0 <= j && j < derived().M && 0 <= k && k < derived().L) {return true;}
                    return ausaxs::detail::report_index_3d("Indexer3D", i, j, k, derived().N, derived().M, derived().L);
                }() && "Indexer3D: Index out of bounds.");
                return derived().data[k + derived().L * (j + derived().M * i)]; 
            }

            constexpr const auto& linear_index(int i) const { 
                assert([&]() -> bool {
                    if (0 <= i && i < derived().N*derived().M*derived().L) {return true;}
                    return ausaxs::detail::report_index_1d("Indexer3D::linear_index", i, derived().N*derived().M*derived().L);
                }() && "Indexer3D::linear_index: Index out of bounds.");
                return derived().data[i]; 
            }

            constexpr const auto& linear_index(int ij, int k) const { 
                return linear_index(ij * derived().L + k);
            }

            constexpr auto& linear_index(int i) { 
                assert([&]() -> bool {
                    if (0 <= i && i < derived().N*derived().M*derived().L) {return true;}
                    return ausaxs::detail::report_index_1d("Indexer3D::linear_index", i, derived().N*derived().M*derived().L);
                }() && "Indexer3D::linear_index: Index out of bounds.");
                return derived().data[i]; 
            }

            constexpr auto& linear_index(int ij, int k) { 
                return linear_index(ij * derived().L + k);
            }

        private:
            Derived& derived() { return static_cast<Derived&>(*this); }
            const Derived& derived() const { return static_cast<const Derived&>(*this); }
    };
}