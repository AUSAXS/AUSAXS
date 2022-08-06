#pragma once

#include <math/Slice.h>

/**
 * @brief \class ConstSlice.
 * 
 * This class represents read-only slices.
 */
template<typename T>
class ConstSlice : public Slice<T> {
	public:
		/**
		 * @brief Construct a new const Slice. 
		 * 
		 * @param data The data backing the original Matrix.
		 * @param N The number of rows of the original Matrix.
		 * @param M The number of columns of the original Matrix.
		 * @param start The start location in the raw data array.
		 * @param step The number of elements to skip between each index in the raw data array. //? Equal to M?
		 * @param length The total number of elements this Slice can access. //? Equal to N?
		 */
		ConstSlice(const std::vector<T>& data, unsigned int N, unsigned int M, unsigned int start, unsigned int step, unsigned int length) : Slice<T>(N, M, start, step, length), data(data) {}

		/**
		 * @brief Destructor.
		 */
		virtual ~ConstSlice() = default;

		/**
		 * @brief Const indexer in this Slice.
		 * 		  Complexity: O(1)
		 */
		const T& operator[](unsigned int i) const override {
			#if (SAFE_MATH)
                // if (__builtin_expect(i >= this->length, false)) {raise(SIGSEGV);}
                if (__builtin_expect(i >= this->length, false)) {throw std::out_of_range("Error in ConstSlice::operator[]: Index out of range.");}
            #endif
			return data[this->start + i*this->step];
		}

	private: 
		const std::vector<T>& data;

		/**
		 * @brief Random-access iterator implementation for const Slices.
		 */
		class ConstIterator {
			using iterator_category = std::random_access_iterator_tag;
			using difference_type = std::ptrdiff_t;
			using value_type = T const;
			using pointer = T const*;
			using reference = T const&;

			public: 
				ConstIterator() : m_ptr(nullptr), step(0) {}
				ConstIterator(pointer ptr, unsigned int step) : m_ptr(ptr), step(step) {}

				// De-reference operators.
				reference operator*() {return *m_ptr;}
				pointer operator->() {return m_ptr;}

				// Increment/decrement operators.
				ConstIterator& operator++() {m_ptr += step; return *this;}
				ConstIterator& operator--() {m_ptr -= step; return *this;}
				ConstIterator operator++(int) {ConstIterator tmp(*this); m_ptr += step; return tmp;}
				ConstIterator operator--(int) {ConstIterator tmp(*this); m_ptr -= step; return tmp;}

				// Arithmetic operators.
				ConstIterator& operator+=(int n) {m_ptr += n*step; return *this;}
				ConstIterator& operator-=(int n) {m_ptr -= n*step; return *this;}
				ConstIterator operator+(int n) {ConstIterator tmp(m_ptr); tmp += n*step; return tmp;}
				ConstIterator operator-(int n) {ConstIterator tmp(m_ptr); tmp -= n*step; return tmp;}
				int operator-(const ConstIterator& other) {return m_ptr - other.m_ptr;}

				// Comparison operators.
				bool operator==(const ConstIterator& other) const {return m_ptr == other.m_ptr;}
				bool operator!=(const ConstIterator& other) const {return m_ptr != other.m_ptr;}
				bool operator<(const ConstIterator& other) const {return m_ptr < other.m_ptr;}
				bool operator>(const ConstIterator& other) const {return m_ptr > other.m_ptr;}
				bool operator<=(const ConstIterator& other) const {return m_ptr <= other.m_ptr;}
				bool operator>=(const ConstIterator& other) const {return m_ptr >= other.m_ptr;}

				// Swap two iterators.
				void swap(ConstIterator& other) {std::swap(m_ptr, other.m_ptr);}

			private: 
				pointer m_ptr;
				unsigned int step;
		};

	public:
		ConstIterator begin() {return ConstIterator(this->data.data() + this->start, this->step);}
		ConstIterator end() {return ConstIterator(this->data.data() + this->start + this->length*this->step, this->step);}
};

template<typename T>
class ConstRow : public ConstSlice<T> {
	public: 
		ConstRow(const std::vector<T>& data, unsigned int N, unsigned int M, unsigned int row) : ConstSlice<T>(data, N, 1, row*M, 1, M) {}

		ConstRow(const ConstSlice<T> s) : ConstSlice<T>(std::move(s)) {}
};

template<typename T>
class ConstColumn : public ConstSlice<T> {
	public: 
		ConstColumn(const std::vector<T>& data, unsigned int N, unsigned int M, unsigned int col) : ConstSlice<T>(data, 1, M, col, M, N) {}

		ConstColumn(const ConstSlice<T> s) : ConstSlice<T>(std::move(s)) {}
};