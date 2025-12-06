#pragma once

#include <utility/indexers/IndexerConcepts.h>

#if (SAFE_MATH) 
    #include <stdexcept>
    #include <string>
#endif

namespace ausaxs::utility::indexer {
    template<typename Derived>
    class Indexer3D {
        protected:
            constexpr const auto& index(int i, int j, int k) const {
                #if (SAFE_MATH)
                    int N = static_cast<int>(derived().N);
                    int M = static_cast<int>(derived().M);
                    int L = static_cast<int>(derived().L);
                    if (i < 0 || N <= i || j < 0 || M <= j || k < 0 || L <= k) {
                        throw std::out_of_range(
                            "Indexer3D: Index out of bounds "
                            "(" + std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) + ") "
                            "should be less than (" + std::to_string(N) + ", " + std::to_string(M) + ", " + std::to_string(L) + ")"
                        );
                    }
                #endif
                return derived().data[k + derived().L * (j + derived().M * i)]; 
            }

            constexpr auto& index(int i, int j, int k) {
                #if (SAFE_MATH)
                    int N = static_cast<int>(derived().N);
                    int M = static_cast<int>(derived().M);
                    int L = static_cast<int>(derived().L);
                    if (i < 0 || N <= i || j < 0 || M <= j || k < 0 || L <= k) {
                        throw std::out_of_range(
                            "Indexer3D: Index out of bounds "
                            "(" + std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) + ") "
                            "should be less than (" + std::to_string(N) + ", " + std::to_string(M) + ", " + std::to_string(L) + ")"
                        );
                    }
                #endif
                return derived().data[k + derived().L * (j + derived().M * i)]; 
            }

            constexpr const auto& linear_index(int i) const { 
                #if (SAFE_MATH)
                    int size = static_cast<int>(derived().N * derived().M * derived().L);
                    if (i < 0 || size <= i) {
                        throw std::out_of_range(
                            "Indexer3D::linear_index: Index out of bounds "
                            "(" + std::to_string(i) + " should be less than " + std::to_string(size) + ")"
                        );
                    }
                #endif
                return derived().data[i]; 
            }

            constexpr auto& linear_index(int i) { 
                #if (SAFE_MATH)
                    int size = static_cast<int>(derived().N * derived().M * derived().L);
                    if (i < 0 || size <= i) {
                        throw std::out_of_range(
                            "Indexer3D::linear_index: Index out of bounds "
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