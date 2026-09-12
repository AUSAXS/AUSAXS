// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <cassert>

#ifndef NDEBUG
    #include <iostream>  // only the asserts below print
#endif

namespace ausaxs::utility::indexer {
    /**
     * @brief CRTP mixin providing element access for a two-dimensional container.
     *        The deriving class must expose a contiguous @c data member (row-major, with row length @c M) and the dimensions @c N and @c M. 
     */
    template<typename Derived>
    class Indexer2D {
        Indexer2D() = default;
        friend Derived;

        protected:
            constexpr const auto& index(int i, int j) const {
                assert([&]() -> bool {
                    if (0 <= i && i < derived().N && 0 <= j && j < derived().M) {return true;}
                    std::cout << "Indexer2D: Index out of bounds (" << i << ", " << j << ") should be less than (" << derived().N << ", " << derived().M << ")" << std::endl;
                    return false;
                }() && "Indexer2D: Index out of bounds.");
                return derived().data[j + derived().M * i]; 
            }

            constexpr auto& index(int i, int j) {
                assert([&]() -> bool {
                    if (0 <= i && i < derived().N && 0 <= j && j < derived().M) {return true;}
                    std::cout << "Indexer2D: Index out of bounds (" << i << ", " << j << ") should be less than (" << derived().N << ", " << derived().M << ")" << std::endl;
                    return false;
                }() && "Indexer2D: Index out of bounds.");
                return derived().data[j + derived().M * i]; 
            }

            constexpr const auto& linear_index(int i) const { 
                assert([&]() -> bool {
                    if (0 <= i && i < derived().N*derived().M) {return true;}
                    std::cout << "Indexer2D::linear_index: Index out of bounds (" << i << " should be less than " << derived().N*derived().M << ")" << std::endl;
                    return false;
                }() && "Indexer2D::linear_index: Index out of bounds.");
                return derived().data[i];
            }

            constexpr auto& linear_index(int i) { 
                assert([&]() -> bool {
                    if (0 <= i && i < derived().N*derived().M) {return true;}
                    std::cout << "Indexer2D::linear_index: Index out of bounds (" << i << " should be less than " << derived().N*derived().M << ")" << std::endl;
                    return false;
                }() && "Indexer2D::linear_index: Index out of bounds.");
                return derived().data[i];
            }

        private:
            Derived& derived() { return static_cast<Derived&>(*this); }
            const Derived& derived() const { return static_cast<const Derived&>(*this); }
    };
}