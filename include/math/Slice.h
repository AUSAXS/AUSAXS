#pragma once

template<class T>
class Vector;

#include <vector>
#include <utility/Exceptions.h>

#define SAFE_MATH true

/**
 * @brief \class Slice
 * 
 * A Slice is a type of cut of a Matrix, typically either in the form of a single row or column. 
 * The idea with this class is to give mutable access to such a structure. 
 */
template<typename T>
class Slice {
    public:
        /**
         * @brief Constructor. 
         * 
         * @param N The number of rows of this Slice. 
         * @param M The number of columns of this Slice. 
         * @param start The start location in the raw data array. 
         * @param step The number of elements to skip between each index in the raw data array. 
         * @param length The total number of elements this Slice can access. 
         */
        Slice(unsigned int N, unsigned int M, unsigned int start, unsigned int step, unsigned int length) : N(N), M(M), start(start), step(step), length(length) {}

        /**
         * @brief Get the dot product with another Slice. 
		 * 		  Complexity: O(N)
         */
        template<typename Q>
        double dot(const Slice<Q>& s) {return dot(s.operator Vector<T>());}

        /**
         * @brief Get the dot product with a Vector.
		 * 		  Complexity: O(N)
         */
        template<typename Q>
        double dot(const Vector<Q>& v) {return operator Vector<T>().dot(v);}

        /**
         * @brief Get the norm of this Slice. 
		 * 		  Complexity: O(N)
		 */
        double norm() {return operator Vector<T>().norm();}

		/**
		 * @brief Get the length of this Slice.
		 */
		size_t size() const noexcept {return length;}

        /**
         * @brief Cast this Slice into a Vector. 
		 * 		  Complexity: O(N)
         */
        operator Vector<T>() const {
            Vector<T> v(length);
            for (unsigned int i = 0; i < length; i++) {
                v[i] = operator[](i);
            }
            return v;
        }

		/**
		 * @brief Cast this Slice into a std::vector.
		 * 		  Complexity: O(N)
		 */
		operator std::vector<T>() const {
			std::vector<T> v(length);
			for (unsigned int i = 0; i < length; i++) {
				v[i] = operator[](i);
			}
			return v;
		}

		/**
		 * @brief Cast this Slice into a std::vector.
		 * 		  Complexity: O(N)
		 */
		std::vector<T> to_vector() const {
			return operator std::vector<T>();
		}

        /**
         * @brief Mutable indexer in this Slice.
		 * 		  Complexity: O(1)
		 */
        virtual const T& operator[](unsigned int j) const = 0;

		/**
		 * @brief Get the final element in this Slice.
		 */
		const T& back() const {return operator[](length-1);}

        /**
         * @brief Get a string representation of this Slice. 
         */
        std::string to_string() const {return operator Vector<T>().to_string();}

        /**
         * @brief Output the string representation of this vector to a stream. 
         */
        friend std::ostream& operator<<(std::ostream& os, const Slice<T>& s) {os << s.to_string(); return os;}

        const unsigned int N, M;
        unsigned int start, step, length;
};

// Add a Vector to this Slice. 
template<typename T, typename Q>
Vector<T> operator+(const Slice<T>& s, const Vector<Q>& v) {return s.operator Vector<T>().operator+(v);}

// Subtract a Vector from this Slice.
template<typename T, typename Q>
Vector<T> operator-(const Slice<T>& s, const Vector<Q>& v) {return s.operator Vector<T>().operator-(v);}

// Multiply each element of this Slice with a constant. 
template<typename T>
Vector<T> operator*(const Slice<T>& s, const double a) {return s.operator Vector<T>().operator*(a);}

// Divide each element of this Slice with a constant.
template<typename T>
Vector<T> operator/(const Slice<T>& s, const double a) {return s.operator Vector<T>().operator/(a);}

// Check if this Slice is equal to a given Vector.
template<typename T, typename Q>
bool operator==(const Slice<T>& s, const Vector<Q>& v) {return s.operator Vector<T>().operator==(v);}

// Check if this Slice is not equal to a given Vector.
template<typename T, typename Q>
bool operator!=(const Slice<T>& s, const Vector<Q>& v) {return s.operator Vector<T>().operator!=(v);}

// VectorSlice class - unused. 
template<typename T>
class VectorSlice : public Slice<T> {
	public: 
		VectorSlice(const std::vector<T>& data) : Slice<T>(data.size(), 1, 0, 1, data.size()), data(data) {}
		virtual ~VectorSlice() {}

		T& operator[](int i) override;
		const std::vector<T>& data;
};

/**
 * @brief A representation of a mutable slice of a Matrix. 
 */
template<typename T>
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
		T& back() {return this->operator[](this->length-1);}

		/**
		 * @brief Mutable indexer in this Slice.
		 * 		  Complexity: O(1)
		 */
		T& operator[](unsigned int i) {
			#if (SAFE_MATH)
                if (__builtin_expect(i >= this->length, false)) {throw std::out_of_range("Error in Matrix::index: Row index out of range.");}
            #endif
			return data[this->start + i*this->step];
		}

		/**
		 * @brief Const indexer in this Slice.
		 * 		  Complexity: O(1)
		 */
		const T& operator[](unsigned int i) const override {
			#if (SAFE_MATH)
                if (__builtin_expect(i >= this->length, false)) {throw std::out_of_range("Error in Matrix::index: Row index out of range.");}
            #endif
			return data[this->start + i*this->step];
		}

		/**
		 * @brief Assign a Vector to this Slice.
		 * 		  Complexity: O(N)
		 */
		MutableSlice& operator=(const Vector<T>& v) {
			#if (SAFE_MATH)
				if (v.size() != this->length) {
					throw except::invalid_argument("Error in MutableSlice::operator=: Vector of size \"" + std::to_string(v.size()) + "\" does not fit in slice of size \"" + std::to_string(this->length) + "\".");
				}
			#endif

			for (unsigned int i = 0; i < this->length; i++) {
				operator[](i) = v[i];
			}
			return *this;
		}

		/**
		 * @brief Subtract another Slice from this Slice.
		 * 		  Complexity: O(N)
		 */
		template<typename Q>
		MutableSlice& operator-=(const Slice<Q>& s) {return operator-=(s.operator Vector<Q>());}

		/**
		 * @brief Subtract a Vector from this Slice.
		 * 		  Complexity: O(N)
		 */
		template<typename Q>
		MutableSlice& operator-=(const Vector<Q>& v) {
			for (unsigned int i = 0; i < this->length; i++) {
				operator[](i) -= v[i];
			}
			return *this;
		}

		/**
		 * @brief Add another Slice to this Slice.
		 * 		  Complexity: O(N)
		 */
		template<typename Q>
		MutableSlice& operator+=(const Slice<Q>& s) {return operator+=(s.operator Vector<Q>());}

		/**
		 * @brief Add a Vector to this Slice.
		 * 		  Complexity: O(N)
		 */
		template<typename Q>
		MutableSlice& operator+=(const Vector<Q>& v) {
			for (unsigned int i = 0; i < this->length; i++) {
				operator[](i) += v[i];
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
                if (__builtin_expect(i >= this->length, false)) {throw std::out_of_range("Error in Matrix::index: Row index out of range.");}
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

template<typename T>
class Row : public MutableSlice<T> {
	public: 
		Row(std::vector<T>& data, unsigned int N, unsigned int M, unsigned int row) : MutableSlice<T>(data, N, 1, row*M, 1, M) {}

		Row(MutableSlice<T> s) : MutableSlice<T>(std::move(s)) {}

		Row<T>& operator=(const Vector<T>& v) {MutableSlice<T> s(*this); MutableSlice<T>::operator=(v); return *this;}

		Row<T>& operator=(const Row<T>& s) {
			for (unsigned int i = 0; i < this->length; i++) {
				this->operator[](i) = s[i];
			}
			return *this;
		}

		Row<T>& operator=(const ConstRow<T>& s) {
			for (unsigned int i = 0; i < this->length; i++) {
				this->operator[](i) = s[i];
			}
			return *this;
		}
};

template<typename T>
class Column : public MutableSlice<T> {
	public: 
		Column(std::vector<T>& data, unsigned int N, unsigned int M, unsigned int col) : MutableSlice<T>(data, 1, M, col, M, N) {}

		Column(MutableSlice<T> s) : MutableSlice<T>(std::move(s)) {}

		Column<T>& operator=(const Vector<T>& v) {MutableSlice<T> s(*this); MutableSlice<T>::operator=(v); return *this;}

		Column<T>& operator=(const Column<T>& s) {
			for (unsigned int i = 0; i < this->length; i++) {
				this->operator[](i) = s[i];
			}
			return *this;
		}

		Column<T>& operator=(const ConstColumn<T>& s) {
			for (unsigned int i = 0; i < this->length; i++) {
				this->operator[](i) = s[i];
			}
			return *this;
		}
};