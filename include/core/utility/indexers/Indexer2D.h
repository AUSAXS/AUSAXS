#pragma once

#include <utility/indexers/IndexerConcepts.h>

#if (SAFE_MATH) 
    #include <stdexcept>
    #include <string>
#endif

namespace ausaxs::utility::indexer {
    template<ValidIndexable2D Derived>
    class Indexer2D {
        using value_type = Derived::value_type;
        protected:
            const value_type& index_impl(int i, int j) const {
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

            value_type& index_impl(int i, int j) {
                return const_cast<value_type&>(static_cast<const Indexer2D&>(*this).index_impl(i, j));
            }

            const auto& linear_index_impl(int i) const { 
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

            value_type& linear_index_impl(int i) { 
                return const_cast<value_type&>(static_cast<const Indexer2D&>(*this).linear_index_impl(i));
            }

        private:
            Derived& derived() { return static_cast<Derived&>(*this); }
            const Derived& derived() const { return static_cast<const Derived&>(*this); }
    };
}