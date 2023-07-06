#pragma once

#include <concepts>

template <typename T>
class SliceIterator {
    using iterator_category = std::random_access_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = T;
    using pointer = T*;
    using reference = T&;

    public:
        SliceIterator() : m_ptr(nullptr), step(0) {}
        SliceIterator(pointer ptr, unsigned int step) : m_ptr(ptr), step(step) {}

        // Dereference operators.
        reference operator*() { return *m_ptr; }
        pointer operator->() { return m_ptr; }

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
        void swap(SliceIterator& other) { std::swap(m_ptr, other.m_ptr); }

    private:
        pointer m_ptr;
        unsigned int step;
};