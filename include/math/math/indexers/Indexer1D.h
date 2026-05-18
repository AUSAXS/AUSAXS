// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#if (SAFE_MATH) 
    #include <stdexcept>
    #include <string>
#endif

namespace ausaxs::utility::indexer {
    /**
     * @brief CRTP mixin providing element access for a one-dimensional container.
     *
     * The deriving class must expose a contiguous @c data member and a @c size() method. When the
     * SAFE_MATH macro is set, every access is bounds-checked and throws std::out_of_range on failure;
     * otherwise the checks compile away to a plain indexed access.
     */
    template<typename Derived>
    class Indexer1D {
        protected:
            constexpr auto& index(int i) {
                #if (SAFE_MATH)
                    int N = static_cast<int>(derived().size());
                    if (i < 0 || N <= i) {
                        throw std::out_of_range(
                            "Indexer1D: Index out of bounds "
                            "(" + std::to_string(i) + " should be less than " + std::to_string(N) + ")"
                        );
                    }
                #endif
                return derived().data[i]; 
            }

            constexpr const auto& index(int i) const {
                #if (SAFE_MATH)
                    int N = static_cast<int>(derived().size());
                    if (i < 0 || N <= i) {
                        throw std::out_of_range(
                            "Indexer1D: Index out of bounds "
                            "(" + std::to_string(i) + " should be less than " + std::to_string(N) + ")"
                        );
                    }
                #endif
                return derived().data[i];
            }

            constexpr auto& linear_index(int i) { return index(i); }
            constexpr const auto& linear_index(int i) const { return index(i); }

        private:
            Derived& derived() { return static_cast<Derived&>(*this); }
            const Derived& derived() const { return static_cast<const Derived&>(*this); }
    };
}