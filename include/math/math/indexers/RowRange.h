// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <cstddef>
#include <iterator>
#include <span>

namespace ausaxs::utility::indexer {
    /**
     * @brief A range over @a count consecutive rows of length @a L in contiguous storage, each viewed as a std::span.
     *        The iterators only refer to the storage, so they stay valid after the range itself is gone.
     */
    template<typename T>
    class RowRange {
        public:
            class iterator {
                public:
                    using value_type = std::span<T>;
                    using difference_type = std::ptrdiff_t;
                    using iterator_concept = std::forward_iterator_tag;

                    iterator() = default;
                    iterator(T* first, std::size_t L, std::size_t row) : first(first), L(L), current(row) {}

                    std::span<T> operator*() const {return {first + current*L, L};}
                    iterator& operator++() {++current; return *this;}
                    iterator operator++(int) {auto tmp = *this; ++current; return tmp;}
                    bool operator==(const iterator& other) const {return current == other.current;}

                private:
                    T* first = nullptr;
                    std::size_t L = 0, current = 0;
            };

            RowRange(T* first, int count, int L) : first(first), count(static_cast<std::size_t>(count)), L(static_cast<std::size_t>(L)) {}

            iterator begin() const {return {first, L, 0};}
            iterator end() const {return {first, L, count};}

        private:
            T* first;
            std::size_t count, L;
    };
}
