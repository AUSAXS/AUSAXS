#pragma once

template<class T>
class Vector;

#include <vector>
#include <stdexcept>
#include <iostream>

using std::vector;

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
        Slice(int N, int M, int start, int step, int length) : N(N), M(M), start(start), step(step), length(length) {}

        /**
         * @brief Get the dot product with another Slice. 
         */
        template<typename Q>
        double dot(const Slice<Q>& s) {return dot(s.operator Vector<T>());}

        /**
         * @brief Get the dot product with a Vector.
         */
        template<typename Q>
        double dot(const Vector<Q>& v) {return operator Vector<T>().dot(v);}

        /**
         * @brief Get the norm of this Slice. 
         */
        double norm() {return operator Vector<T>().norm();}

        /**
         * @brief Cast this Slice into a Vector. 
         */
        operator Vector<T>() const {
            Vector<T> v(length);
            for (int i = 0; i < length; i++) {
                v[i] = operator[](i);
            }
            return v;
        }

        /**
         * @brief Mutable indexer in this Slice.
         * 
         * @param j The index to access.
         */
        virtual const T& operator[](unsigned int j) const = 0;

        /**
         * @brief Get a string representation of this Slice. 
         */
        std::string to_string() const {return operator Vector<T>().to_string();}

        /**
         * @brief Output the string representation of this vector to a stream. 
         */
        template<typename Q>
        friend std::ostream& operator<<(std::ostream& os, const Slice<Q>& s) {os << s.to_string(); return os;}

        const size_t N, M;
        int start, step, length;
};

// Add a Vector to this Slice. 
template<typename T, typename Q>
Vector<T> operator+(const Slice<T>& s, const Vector<Q>& v);

// Subtract a Vector from this Slice.
template<typename T, typename Q>
Vector<T> operator-(const Slice<T>& s, const Vector<Q>& v);

// Multiply each element of this Slice with a constant. 
template<typename T>
Vector<T> operator*(const Slice<T>& s, const double a);

// Divide each element of this Slice with a constant.
template<typename T>
Vector<T> operator/(const Slice<T>& s, const double a);

// Check if this Slice is equal to a given Vector.
template<typename T, typename Q>
bool operator==(const Slice<T>& s, const Vector<Q>& v);

// Check if this Slice is not equal to a given Vector.
template<typename T, typename Q>
bool operator!=(const Slice<T>& s, const Vector<Q>& v);

// VectorSlice class - unused. 
template<typename T>
class VectorSlice : public Slice<T> {
	public: 
		VectorSlice(const vector<T>& data) : Slice<T>(data.size(), 1, 0, 1, data.size()), data(data) {}
		virtual ~VectorSlice() {}

		T& operator[](int i) override;
		const vector<T>& data;
};

/**
 * @brief \class MutableSlice. 
 * 
 * This class represents a mutable slice. 
 */
template<typename T>
class MutableSlice : public Slice<T> {
	public:
		MutableSlice(vector<T>& data, int N, int M, int start, int step, int length) : Slice<T>(N, M, start, step, length), data(data) {}
	
		virtual ~MutableSlice() = default;

		T& operator[](unsigned int i) {return data[this->start + i*this->step];}

		const T& operator[](unsigned int i) const override {return data[this->start + i*this->step];}

		MutableSlice& operator=(const Vector<T>& v) {
			for (int i = 0; i < this->length; i++) {
				operator[](i) = v[i];
			}
			return *this;
		}

		template<typename Q>
		MutableSlice& operator-=(const Slice<Q>& s) {return operator-=(s.operator Vector<Q>());}

		template<typename Q>
		MutableSlice& operator-=(const Vector<Q>& v) {
			for (int i = 0; i < this->length; i++) {
				operator[](i) -= v[i];
			}
			return *this;
		}

		template<typename Q>
		MutableSlice& operator+=(const Slice<Q>& s) {return operator+=(s.operator Vector<Q>());}

		template<typename Q>
		MutableSlice& operator+=(const Vector<Q>& v) {
			for (int i = 0; i < this->length; i++) {
				operator[](i) += v[i];
			}
			return *this;
		}

		vector<T>& data;
};

/**
 * @brief \class ConstSlice.
 * 
 * This class represents read-only slices.
 */
template<typename T>
class ConstSlice : public Slice<T> {
	public:
		ConstSlice(const vector<T>& data, int N, int M, int start, int step, int length) : Slice<T>(N, M, start, step, length), data(data) {}

		virtual ~ConstSlice() = default;

		const T& operator[](unsigned int i) const override {return data[this->start + i*this->step];}

		const vector<T>& data;
};

template<typename T>
class ConstRow : public ConstSlice<T> {
	public: 
		ConstRow(const vector<T>& data, int N, int M, int row) : ConstSlice<T>(data, N, 1, row*M, 1, M) {}

		ConstRow(const ConstSlice<T>& s) : ConstSlice<T>(std::move(s)) {}
};

template<typename T>
class ConstColumn : public ConstSlice<T> {
	public: 
		ConstColumn(const vector<T>& data, int N, int M, int col) : ConstSlice<T>(data, 1, M, col, M, N) {}

		ConstColumn(const ConstSlice<T>& s) : ConstSlice<T>(std::move(s)) {}
};

template<typename T>
class Row : public MutableSlice<T> {
	public: 
		Row(vector<T>& data, int N, int M, int row) : MutableSlice<T>(data, N, 1, row*M, 1, M) {}

		Row(MutableSlice<T>& s) : MutableSlice<T>(std::move(s)) {}

		Row<T>& operator=(const Vector<T>& v) {MutableSlice<T> s(*this); MutableSlice<T>::operator=(v); return *this;}

		Row<T>& operator=(const Row<T>& s) {
			for (int i = 0; i < this->length; i++) {
				this->operator[](i) = s[i];
			}
			return *this;
		}

		Row<T>& operator=(const ConstRow<T>& s) {
			for (int i = 0; i < this->length; i++) {
				this->operator[](i) = s[i];
			}
			return *this;
		}
};

template<typename T>
class Column : public MutableSlice<T> {
	public: 
		Column(vector<T>& data, int N, int M, int col) : MutableSlice<T>(data, 1, M, col, M, N) {}

		Column(MutableSlice<T>& s) : MutableSlice<T>(std::move(s)) {}

		Column<T>& operator=(const Vector<T>& v) {MutableSlice<T> s(*this); MutableSlice<T>::operator=(v); return *this;}

		Column<T>& operator=(const Column<T>& s) {
			for (int i = 0; i < this->length; i++) {
				this->operator[](i) = s[i];
			}
			return *this;
		}

		Column<T>& operator=(const ConstColumn<T>& s) {
			for (int i = 0; i < this->length; i++) {
				this->operator[](i) = s[i];
			}
			return *this;
		}
};

template<typename T, typename Q>
Vector<Q> operator+(const Slice<T>& s, const Vector<Q>& v) {return s.operator Vector<T>().operator+(v);}

template<typename T, typename Q>
Vector<Q> operator-(const Slice<T>& s, const Vector<Q>& v) {return s.operator Vector<T>().operator-(v);}

template<typename T>
Vector<T> operator*(const Slice<T>& s, const double a) {return s.operator Vector<T>().operator*(a);}

template<typename T>
Vector<T> operator/(const Slice<T>& s, const double a) {return s.operator Vector<T>().operator/(a);}

template<typename T, typename Q>
bool operator==(const Slice<T>& s, const Vector<Q>& v) {return s.operator Vector<T>().operator==(v);}

template<typename T, typename Q>
bool operator!=(const Slice<T>& s, const Vector<Q>& v) {return s.operator Vector<T>().operator!=(v);}