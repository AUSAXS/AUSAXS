#pragma once

#if (SAFE_MATH) 
    #include <stdexcept>
    #include <string>
#endif

namespace ausaxs::utility::indexer {
    template<typename Derived>
    class Indexer2D {
        protected:
            constexpr const auto& index(int i, int j) const {
                #if (SAFE_MATH)
                    int N = static_cast<int>(derived().N);
                    int M = static_cast<int>(derived().M);
                    if (i < 0 || N <= i || j < 0 || M <= j) {
                        throw std::out_of_range(
                            "Indexer2D: Index out of bounds "
                            "(" + std::to_string(i) + ", " + std::to_string(j) + ") "
                            "should be less than (" + std::to_string(N) + ", " + std::to_string(M) + ")"
                        );
                    }
                #endif
                return derived().data[j + derived().M * i]; 
            }

            constexpr auto& index(int i, int j) {
                #if (SAFE_MATH)
                    int N = static_cast<int>(derived().N);
                    int M = static_cast<int>(derived().M);
                    if (i < 0 || N <= i || j < 0 || M <= j) {
                        throw std::out_of_range(
                            "Indexer2D: Index out of bounds "
                            "(" + std::to_string(i) + ", " + std::to_string(j) + ") "
                            "should be less than (" + std::to_string(N) + ", " + std::to_string(M) + ")"
                        );
                    }
                #endif
                return derived().data[j + derived().M * i]; 
            }

            constexpr const auto& linear_index(int i) const { 
                #if (SAFE_MATH)
                    int size = static_cast<int>(derived().N * derived().M);
                    if (i < 0 || size <= i) {
                        throw std::out_of_range(
                            "Indexer2D::linear_index: Index out of bounds "
                            "(" + std::to_string(i) + " should be less than " + std::to_string(size) + ")"
                        );
                    }
                #endif
                return derived().data[i];
            }

            constexpr auto& linear_index_impl(int i) { 
                #if (SAFE_MATH)
                    int size = static_cast<int>(derived().N * derived().M);
                    if (i < 0 || size <= i) {
                        throw std::out_of_range(
                            "Indexer2D::linear_index: Index out of bounds "
                            "(" + std::to_string(i) + " should be less than " + std::to_string(size) + ")"
                        );
                    }
                #endif
                return derived().data[i];
            }

        private:
            Derived& derived() { return static_cast<Derived&>(*this); }
            const Derived& derived() const { return static_cast<const Derived&>(*this); }
    };
}