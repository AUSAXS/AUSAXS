#pragma once

#include <math/slices/Slice.h>

#include <concepts>

// define slice concepts
// template<typename T>
// concept ConstSliceConcept = std::derived_from<T, Slice<T>>;

// template<typename T>
// concept MutableSliceConcept = std::derived_from<T, MutableSlice<T>>;

// template<ConstSliceConcept T>
// struct ConstOps {
// 	// double ops
// 	T operator+(double rhs) const {
// 		T result = *this;
// 		for (unsigned int i = 0; i < result.size(); i++) {
// 			result[i] += rhs;
// 		}
// 		return result;
// 	}
// 	T operator-(double rhs) const {return operator+(-rhs);}

// 	T operator*(double rhs) const {
// 		T result = *this;
// 		for (unsigned int i = 0; i < result.size(); i++) {
// 			result[i] *= rhs;
// 		}
// 		return result;
// 	}
// 	T operator/(double rhs) const {return operator*(1/rhs);}

// 	// slice ops
// 	template<ConstSliceConcept Q>
// 	T& operator+=(const Q& rhs) {return *this = *this + rhs;}
// };

// template<MutableSliceConcept T>
// struct MutableOps : public ConstOps<T> {
// 	// double ops
// 	T& operator+=(double rhs) {return *this = *this + rhs;}
// 	T& operator-=(double rhs) {return *this = *this - rhs;}
// 	T& operator*=(double rhs) {return *this = *this * rhs;}
// 	T& operator/=(double rhs) {return *this = *this / rhs;}
// };

/**
 * @brief A mutable slice of a Matrix. 
 */
template<numeric T>
class MutableSlice : public Slice<T> {
	public:
		/**
		 * @brief Construct a new mutable Slice. 
		 * 
		 * @param data The data backing the original Matrix.
		 * @param N The number of rows of the original Matrix.
		 * @param M The number of columns of the original Matrix.
		 * @param start The start location in the raw data array.
		 * @param step The number of elements to skip between each index in the raw data array. //? Equal to M?
		 * @param length The total number of elements this Slice can access. //? Equal to N?
		 */
		MutableSlice(std::vector<T>& data, unsigned int N, unsigned int M, unsigned int start, unsigned int step, unsigned int length) : Slice<T>(N, M, start, step, length), data(data) {}

		/**
		 * @brief Destructor.
		 */
		virtual ~MutableSlice() = default;

		/**
		 * @brief Get the final element in this Slice.
		 */
		T& back() {
			#if SAFE_MATH
				if (this->size() == 0) [[unlikely]] {throw std::out_of_range("MutableSlice::back(): Slice is empty.");}
			#endif

			return (*this)[this->length-1];
		}

		/**
		 * @brief Get the first element in this Slice.
		 */
		T& front() {
			#if SAFE_MATH
				if (this->size() == 0) [[unlikely]] {throw std::out_of_range("MutableSlice::front(): Slice is empty.");}
			#endif

			return (*this)[0];
		}

		T& operator[](unsigned int i) {
			#if (SAFE_MATH)
                if (i >= this->length) [[unlikely]] {throw std::out_of_range("MutableSlice::operator[]: Index out of range.");}
            #endif
			return data[this->start + i*this->step];
		}

		const T& operator[](unsigned int i) const override {
			#if (SAFE_MATH)
                if (i >= this->length) [[unlikely]] {throw std::out_of_range("MutableSlice::operator[]: Index out of range.");}
            #endif
			return data[this->start + i*this->step];
		}

		MutableSlice<T>& operator=(const Vector<T>& v) {
			#if (SAFE_MATH)
				if (v.size() != this->length) [[unlikely]] {
					throw std::invalid_argument("MutableSlice::operator=: Vector of size \"" + std::to_string(v.size()) + "\" does not fit in slice of size \"" + std::to_string(this->length) + "\".");
				}
			#endif

			for (unsigned int i = 0; i < this->length; i++) {
				(*this)[i] = v[i];
			}
			return *this;
		}

		template<typename Q>
		MutableSlice<T>& operator-=(const Slice<Q>& s) {
			for (unsigned int i = 0; i < this->length; i++) {
				(*this)[i] -= s[i];
			}
			return *this;
		}

		template<typename Q>
		MutableSlice<T>& operator-=(const Vector<Q>& v) {
			for (unsigned int i = 0; i < this->length; i++) {
				(*this)[i] -= v[i];
			}
			return *this;
		}

		template<typename Q>
		MutableSlice<T>& operator-=(double a) {
			return operator+=(-a);
		}

		template<typename Q>
		MutableSlice<T>& operator+=(const Slice<Q>& s) {return operator+=(s.operator Vector<Q>());}


		template<typename Q>
		MutableSlice<T>& operator+=(double a) {
			for (unsigned int i = 0; i < this->length; i++) {
				(*this)[i] += a;
			}
			return *this;
		}

		template<typename Q>
		MutableSlice<T>& operator+=(const Vector<Q>& v) {
			for (unsigned int i = 0; i < this->length; i++) {
				(*this)[i] += v[i];
			}
			return *this;
		}

		MutableSlice<T>& operator*=(double a) {
			for (unsigned int i = 0; i < this->length; i++) {
				(*this)[i] *= a;
			}
			return *this;
		}

		MutableSlice<T>& operator/=(double a) {
			for (unsigned int i = 0; i < this->length; i++) {
				(*this)[i] /= a;
			}
			return *this;
		}

	private: 
		std::vector<T>& data;

		/**
		 * @brief Random-access iterator implementation for mutable Slices.
		 */
		class Iterator {
			using iterator_category = std::random_access_iterator_tag;
			using difference_type = std::ptrdiff_t;
			using value_type = T;
			using pointer = T*;
			using reference = T&;

			public: 
				Iterator() : m_ptr(nullptr), step(0) {}
				Iterator(pointer ptr, unsigned int step) : m_ptr(ptr), step(step) {}

				// De-reference operators.
				reference operator*() {return *m_ptr;}
				pointer operator->() {return m_ptr;}

				// Increment/decrement operators.
				Iterator& operator++() {m_ptr += step; return *this;}
				Iterator& operator--() {m_ptr -= step; return *this;}
				Iterator operator++(int) {Iterator tmp(*this); m_ptr += step; return tmp;}
				Iterator operator--(int) {Iterator tmp(*this); m_ptr -= step; return tmp;}

				// Arithmetic operators.
				Iterator& operator+=(int n) {m_ptr += n*step; return *this;}
				Iterator& operator-=(int n) {m_ptr -= n*step; return *this;}
				Iterator operator+(int n) {Iterator tmp(m_ptr); tmp += n*step; return tmp;}
				Iterator operator-(int n) {Iterator tmp(m_ptr); tmp -= n*step; return tmp;}
				int operator-(const Iterator& other) {return m_ptr - other.m_ptr;}

				// Comparison operators.
				bool operator==(const Iterator& other) const {return m_ptr == other.m_ptr;}
				bool operator!=(const Iterator& other) const {return m_ptr != other.m_ptr;}
				bool operator<(const Iterator& other) const {return m_ptr < other.m_ptr;}
				bool operator>(const Iterator& other) const {return m_ptr > other.m_ptr;}
				bool operator<=(const Iterator& other) const {return m_ptr <= other.m_ptr;}
				bool operator>=(const Iterator& other) const {return m_ptr >= other.m_ptr;}

				// Swap two iterators.
				void swap(Iterator& other) {std::swap(m_ptr, other.m_ptr);}

			private: 
				pointer m_ptr;
				unsigned int step;
		};
	
	public: 		
		Iterator begin() {return Iterator(this->data.data() + this->start, this->step);}
		Iterator end() {return Iterator(this->data.data() + this->start + this->length*this->step, this->step);}
};

// Add a Vector to this Slice. 
template<typename T, typename Q>
MutableSlice<T> operator+(MutableSlice<T> left, const Vector<Q>& right) {return left += right;}

// Subtract a Vector from this Slice.
template<typename T, typename Q>
MutableSlice<T> operator-(MutableSlice<T> left, const Vector<Q>& right) {return left -= right;}

// Add a constant to each element of this Slice. 
template<typename T>
MutableSlice<T> operator+(MutableSlice<T> left, double right) {
	for (auto& e : left) {e += right;}
	return left;
}

// Subtract a constant to each element of this Slice. 
template<typename T>
MutableSlice<T> operator-(MutableSlice<T> left, double right) {return left + -right;}

// Multiply each element of this Slice with a constant. 
template<typename T>
MutableSlice<T> operator*(MutableSlice<T> left, double right) {return left *= right;}

// Divide each element of this Slice with a constant.
template<typename T>
MutableSlice<T> operator/(MutableSlice<T> left, double right) {return left /= right;}

// Check if this Slice is equal to a given Vector.
template<typename T, typename Q>
bool operator==(const MutableSlice<T>& s, const Vector<Q>& v) {return s.operator Vector<T>().operator==(v);}

// Check if this Slice is not equal to a given Vector.
template<typename T, typename Q>
bool operator!=(const MutableSlice<T>& s, const Vector<Q>& v) {return !(s == v);}

template<typename T>
class Row : public MutableSlice<T> {
	public: 
		Row(std::vector<T>& data, unsigned int N, unsigned int M, unsigned int row) : MutableSlice<T>(data, N, 1, row*M, 1, M) {}

		Row(MutableSlice<T> s) : MutableSlice<T>(std::move(s)) {}

		// Vector assignment
		Row<T>& operator=(const Vector<T>& v) {
			#if (SAFE_MATH)
                if (v.size() != this->size()) [[unlikely]] {
					throw std::invalid_argument("Column::operator=: Vector of size \"" + std::to_string(v.size()) + "\" does not fit in slice of size \"" + std::to_string(this->size()) + "\".");
				}
            #endif

			MutableSlice<T>::operator=(v); 
			return *this;
		}

		// Row assignment
		Row<T>& operator=(const Row<T>& s) {
			#if (SAFE_MATH)
                if (s.size() != this->size()) [[unlikely]] {
					throw std::invalid_argument("Column::operator=: Slice of size \"" + std::to_string(s.size()) + "\" does not fit in slice of size \"" + std::to_string(this->size()) + "\".");
				}
            #endif

			for (unsigned int i = 0; i < this->length; i++) {
				(*this)[i] = s[i];
			}
			return *this;
		}

		// ConstRow assignment
		template<typename Q, std::enable_if_t<std::is_base_of_v<Slice<T>, Q>, int> = 0>
		Row<T>& operator=(const Q& s) {
			#if (SAFE_MATH)
                if (s.size() != this->size()) [[unlikely]] {
					throw std::invalid_argument("Column::operator=: Slice of size \"" + std::to_string(s.size()) + "\" does not fit in slice of size \"" + std::to_string(this->size()) + "\".");
				}
            #endif

			for (unsigned int i = 0; i < this->length; i++) {
				(*this)[i] = s[i];
			}
			return *this;
		}
};

template<typename T>
class Column : public MutableSlice<T> {
	public: 
		Column(std::vector<T>& data, unsigned int N, unsigned int M, unsigned int col) : MutableSlice<T>(data, 1, M, col, M, N) {}

		Column(MutableSlice<T> s) : MutableSlice<T>(std::move(s)) {}

		// Vector assignment
		Column<T>& operator=(const Vector<T>& v) {
			#if (SAFE_MATH)
                if (v.size() != this->size()) [[unlikely]] {
					throw std::invalid_argument("Column::operator=: Vector of size \"" + std::to_string(v.size()) + "\" does not fit in slice of size \"" + std::to_string(this->size()) + "\".");
				}
            #endif

			MutableSlice<T>::operator=(v);
			return *this;
		}

		// Column assignment
		Column<T>& operator=(const Column<T>& s) {
			#if (SAFE_MATH)
                if (s.size() != this->size()) [[unlikely]] {
					throw std::invalid_argument("Column::operator=: Slice of size \"" + std::to_string(s.size()) + "\" does not fit in slice of size \"" + std::to_string(this->size()) + "\".");
				}
            #endif

			for (unsigned int i = 0; i < this->length; i++) {
				(*this)[i] = s[i];
			}
			return *this;
		}

		// ConstColumn assignment
		template<typename Q, std::enable_if_t<std::is_base_of_v<Slice<T>, Q>, int> = 0>
		Column<T>& operator=(const Q& s) {
			#if (SAFE_MATH)
                if (s.size() != this->size()) [[unlikely]] {
					throw std::invalid_argument("Column::operator=: Slice of size \"" + std::to_string(s.size()) + "\" does not fit in slice of size \"" + std::to_string(this->size()) + "\".");
				}
            #endif

			for (unsigned int i = 0; i < this->length; i++) {
				(*this)[i] = s[i];
			}
			return *this;
		}

		using MutableSlice<T>::operator+=;
		Column<T>& operator+=(double a) {
			for (unsigned int i = 0; i < this->length; i++) {
				(*this)[i] += a;
			}
			return *this;
		}
};