// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <cassert>

#ifndef NDEBUG
    #include <iostream>  // only the asserts below print
#endif

namespace ausaxs::utility::indexer {
    /**
     * @brief CRTP mixin providing element access for a one-dimensional container.
     *        The deriving class must expose a contiguous @c data member and a @c size() method. 
     */
    template<typename Derived>
    class Indexer1D {
        Indexer1D() = default;
        friend Derived;

        protected:
            constexpr auto& index(int i) {
                assert([&]() -> bool {
                    if (0 <= i && i < static_cast<int>(derived().size())) {return true;}
                    std::cout << "Indexer1D: Index out of bounds (" << i << " should be less than " << static_cast<int>(derived().size()) << ")" << std::endl;
                    return false;
                }() && "Indexer1D: Index out of bounds.");
                return derived().data[i]; 
            }

            constexpr const auto& index(int i) const {
                assert([&]() -> bool {
                    if (0 <= i && i < static_cast<int>(derived().size())) {return true;}
                    std::cout << "Indexer1D: Index out of bounds (" << i << " should be less than " << static_cast<int>(derived().size()) << ")" << std::endl;
                    return false;
                }() && "Indexer1D: Index out of bounds.");
                return derived().data[i];
            }

            constexpr auto& linear_index(int i) { return index(i); }
            constexpr const auto& linear_index(int i) const { return index(i); }

        private:
            Derived& derived() { return static_cast<Derived&>(*this); }
            const Derived& derived() const { return static_cast<const Derived&>(*this); }
    };
}