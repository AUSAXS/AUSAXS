#pragma once

#include <utility/indexers/IndexerConcepts.h>

#if (SAFE_MATH) 
    #include <stdexcept>
    #include <string>
#endif

namespace ausaxs::utility::indexer {
    template<ValidIndexable1D Derived>
    class Indexer1D {
        using value_type = Derived::value_type;
        protected:
            const value_type& index_impl(int i) const {
                #if (SAFE_MATH)
                    int N = static_cast<int>(derived().N);
                    if (i < 0 || N <= i) {
                        throw std::out_of_range(
                            "Indexer1D: Index out of bounds "
                            "(" + std::to_string(i) + " should be less than " + std::to_string(N) + ")"
                        );
                    }
                #endif
                return derived().data[i]; 
            }

            value_type& index_impl(int i) {
                return const_cast<value_type&>(static_cast<const Indexer1D&>(*this).index_impl(i));
            }

            value_type linear_index_impl(int i) { return index_impl(i); }
            const value_type& linear_index_impl(int i) const { return index_impl(i); }

        private:
            Derived& derived() { return static_cast<Derived&>(*this); }
            const Derived& derived() const { return static_cast<const Derived&>(*this); }
    };
}