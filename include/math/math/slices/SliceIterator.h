// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <cstddef>
#include <iterator>

namespace ausaxs {
    template <typename T>
    class SliceIterator {
        public:
            using iterator_category = std::forward_iterator_tag;
            using difference_type = std::ptrdiff_t;
            using value_type = T;
            using pointer = T*;
            using reference = T&;

            SliceIterator() : m_ptr(nullptr), step(0) {}
            SliceIterator(pointer ptr, int step) : m_ptr(ptr), step(step) {}

            reference operator*() const { return *m_ptr; }
            pointer operator->() const { return m_ptr; }

            // Increment/decrement operators.
            SliceIterator& operator++() { m_ptr += step; return *this; }
            SliceIterator& operator--() { m_ptr -= step; return *this; }
            SliceIterator operator++(int) { SliceIterator tmp(*this); m_ptr += step; return tmp; }
            SliceIterator operator--(int) { SliceIterator tmp(*this); m_ptr -= step; return tmp; }

            // Arithmetic operators.
            SliceIterator& operator+=(int n) { m_ptr += n * step; return *this; }
            SliceIterator& operator-=(int n) { m_ptr -= n * step; return *this; }
            SliceIterator operator+(int n) { SliceIterator tmp(m_ptr); tmp += n * step; return tmp; }
            SliceIterator operator-(int n) { SliceIterator tmp(m_ptr); tmp -= n * step; return tmp; }
            int operator-(const SliceIterator& other) { return m_ptr - other.m_ptr; }

            // Comparison operators.
            bool operator==(const SliceIterator& other) const { return m_ptr == other.m_ptr; }
            bool operator!=(const SliceIterator& other) const { return m_ptr != other.m_ptr; }
            bool operator<(const SliceIterator& other) const { return m_ptr < other.m_ptr; }
            bool operator>(const SliceIterator& other) const { return m_ptr > other.m_ptr; }
            bool operator<=(const SliceIterator& other) const { return m_ptr <= other.m_ptr; }
            bool operator>=(const SliceIterator& other) const { return m_ptr >= other.m_ptr; }

            // Swap two iterators.
            void swap(SliceIterator& other) noexcept { std::swap(m_ptr, other.m_ptr); }

        private:
            pointer m_ptr;
            int step;
    };
}