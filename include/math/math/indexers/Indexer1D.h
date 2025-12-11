// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#if (SAFE_MATH) 
    #include <stdexcept>
    #include <string>
#endif

namespace ausaxs::utility::indexer {
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