#pragma once

#include <Symbols.h>
#include <utility/Concepts.h>
#include <math/slices/SliceIterator.h>
#include <math/Vector.h>

#include <vector>
#include <initializer_list>

template<numeric T, container_type Container> 
class Slice {
    public: 
        /**
         * @param data The raw data array.
         * @param offset The offset in the raw data array. 
         * @param step The number of elements to skip between each index in the raw data array. 
         * @param length The total number of elements this Slice can access. 
         */
        Slice(Container data, unsigned int offset, unsigned int step, unsigned int length) : data(data), offset(offset), step(step), length(length) {}
        Slice(Slice&& s) noexcept : data(std::move(s.data)), offset(s.offset), step(s.step), length(s.length) {}
        Slice(const Slice& s) : data(s.data), offset(s.offset), step(s.step), length(s.length) {}

        unsigned int size() const noexcept {return length;}

        /**
         * @brief Immutable indexer in this Slice.
		 * 		  Complexity: O(1)
		 */
        virtual const T& operator[](unsigned int j) const {
            #if SAFE_MATH
                if (j >= size()) {throw std::out_of_range("Slice::operator[]: Index out of range.");}
            #endif
            return data[offset + j*step];
        }

		/**
		 * @brief Get the final element in this Slice.
		 */
		const T& back() const {
			#if SAFE_MATH
				if (size() == 0) {throw std::out_of_range("Slice::back(): Slice is empty.");}
			#endif
			return (*this)[length-1];
		}

		/**
		 * @brief Get the final element in this Slice.
		 */
		const T& first() const {
			#if SAFE_MATH
				if (size() == 0) {throw std::out_of_range("Slice::first(): Slice is empty.");}
			#endif
			return (*this)[0];
		}

        template<container_type Q>
        double dot(const Q& rhs) const {
            validate_sizes(rhs.size());
            double sum = 0;
            for (unsigned int i = 0; i < size(); i++) {
                sum += (*this)[i] * rhs[i];
            }
            return sum;
        }

        double norm() {
            return sqrt(dot(*this));
        }

        template<container_type Q>
        bool operator==(const Q& rhs) const {
            validate_sizes(rhs.size());
            bool equal = true;
            for (unsigned int i = 0; i < size(); i++) {
                equal = equal && (*this)[i] == rhs[i];
            }
            return equal;
        }

        template<typename Q, container_type R>
        std::vector<T> operator+(const Slice<Q, R>& rhs) const {
            std::vector<T> result(size());
            for (unsigned int i = 0; i < size(); i++) {
                result[i] = (*this)[i] + rhs[i];
            }
            return result;
        }

        template<typename Q, container_type R>
        std::vector<T> operator-(const Slice<Q, R>& rhs) const {
            std::vector<T> result(size());
            for (unsigned int i = 0; i < size(); i++) {
                result[i] = (*this)[i] - rhs[i];
            }
            return result;
        }

        std::vector<T> operator/(double rhs) const {
            std::vector<T> result(size());
            for (unsigned int i = 0; i < size(); i++) {
                result[i] = (*this)[i] / rhs;
            }
            return result;
        }

        std::vector<T> operator*(double rhs) const {
            std::vector<T> result(size());
            for (unsigned int i = 0; i < size(); i++) {
                result[i] = (*this)[i] * rhs;
            }
            return result;
        }

        operator std::vector<T>() const {
            std::vector<T> v(size());
            for (unsigned int i = 0; i < size(); i++) {
                v[i] = (*this)[i];
            }
            return v;
        }

        operator Vector<T>() const {
            Vector<T> v(size());
            for (unsigned int i = 0; i < size(); i++) {
                v[i] = (*this)[i];
            }
            return v;
        }

        SliceIterator<const T> begin() const {return SliceIterator<const T>(&data[offset], step);}
        SliceIterator<const T> end() const {return SliceIterator<const T>(&data[offset + length*step], step);}

    protected:
        Container data;         // the raw data array backing this Slice
        unsigned int offset;    // offset index in linear data array 
        unsigned int step;      // step size between each index in linear data array
        unsigned int length;    // total number of elements in this Slice

        void validate_sizes(unsigned int length) const {
            #if SAFE_MATH
                if (size() != length) {
                    throw std::invalid_argument("Slice::validate_sizes: Slice of size \"" + std::to_string(length) + "\" does not fit in slice of size \"" + std::to_string(size()) + "\".");
                }
            #endif
        }
};

template<numeric T>
class ConstSlice : public Slice<T, const std::vector<T>&> {
        using Slice<T, const std::vector<T>&>::Slice;
};

template<numeric T>
class MutableSlice : public Slice<T, std::vector<T>&> {
    using data_type = std::vector<T>&;
    public:
        using Slice<T, data_type>::Slice;

        operator ConstSlice<T>() const {
            return ConstSlice<T>(this->data, this->offset, this->step, this->length);
        }

        /**
         * @brief Mutable indexer in this Slice.
         *           Complexity: O(1)
         */
        virtual T& operator[](unsigned int j) {
            #if SAFE_MATH
                if (j >= this->size()) {throw std::out_of_range("Slice::operator[]: Index out of range.");}
            #endif
            return this->data[this->offset + j*this->step];
        }
        using Slice<T, data_type>::operator[];

        MutableSlice& operator=(std::initializer_list<T> rhs) {
            std::vector<T> v(rhs);
            *this = v;
            return *this;
        }

        template<container_type Q>
        MutableSlice& operator=(const Q& rhs) {
            this->validate_sizes(rhs.size());
            for (unsigned int i = 0; i < this->size(); i++) {
                (*this)[i] = rhs[i];
            }
            return *this;
        }

        template<numeric Q, container_type R>
        MutableSlice& operator+=(const Slice<Q, R>& rhs) {
            this->validate_sizes(rhs.size());
            for (unsigned int i = 0; i < this->size(); i++) {
                (*this)[i] += rhs[i];
            }
            return *this;
        }

        template<numeric Q, container_type R>
        MutableSlice& operator-=(const Slice<Q, R>& rhs) {
            this->validate_sizes(rhs.size());
            for (unsigned int i = 0; i < this->size(); i++) {
                (*this)[i] -= rhs[i];
            }
            return *this;
        }


		/**
		 * @brief Get the final element in this Slice.
		 */
		T& back() {
			#if SAFE_MATH
				if (this->size() == 0) {throw std::out_of_range("MutableSlice::back(): Slice is empty.");}
			#endif
			return (*this)[this->size()-1];
		}
        using Slice<T, data_type>::back;

		/**
		 * @brief Get the final element in this Slice.
		 */
		T& first() {
			#if SAFE_MATH
				if (this->size() == 0) {throw std::out_of_range("MutableSlice::first(): Slice is empty.");}
			#endif
			return (*this)[0];
		}
        using Slice<T, data_type>::first;

        SliceIterator<T> begin() {return SliceIterator<T>(&this->data[this->offset], this->step);}
        SliceIterator<T> end() {return SliceIterator<T>(&this->data[this->offset + this->length*this->step], this->step);}
};

template<numeric T>
class ConstRow : public ConstSlice<T> {
    public:
        /**
         * @brief Constructor. 
         * 
         * @param data The raw data array.
         * @param N Unused parameter (for consistency with other constructors).
         * @param M The number of columns of this Slice. 
         * @param row The row index of this ConstRow.
         */
        ConstRow(const std::vector<T>& data, unsigned int, unsigned int M, unsigned int row) : ConstSlice<T>(data, row*M, 1, M) {}
        using ConstSlice<T>::ConstSlice;
};

template<numeric T>
class ConstColumn : public ConstSlice<T> {
    public:
        /**
         * @brief Constructor. 
         * 
         * @param data The raw data array.
         * @param offset The offset in the raw data array. 
         * @param step The number of elements to skip between each index in the raw data array. 
         * @param length The total number of elements this Slice can access. 
         */
        ConstColumn(const std::vector<T>& data, unsigned int N, unsigned int M, unsigned int col) : ConstSlice<T>(data, col, M, N) {}
        using ConstSlice<T>::ConstSlice;
};

template<numeric T>
class MutableRow : public MutableSlice<T> {
    public:
        /**
         * @brief Constructor. 
         * 
         * @param data The raw data array.
         * @param N Unused parameter (for consistency with other constructors).
         * @param M The number of columns of this Slice. 
         * @param row The row index of this ConstRow.
         */
        MutableRow(std::vector<T>& data, unsigned int, unsigned int M, unsigned int row) : MutableSlice<T>(data, row*M, 1, M) {}
        using MutableSlice<T>::MutableSlice;

        // We have to explicitly define this to avoid ambiguity
        MutableRow& operator=(const MutableRow& rhs) {MutableSlice<T>::operator=(rhs); return *this;}
        using MutableSlice<T>::operator=;
};

template<numeric T>
class MutableColumn : public MutableSlice<T> {
    public:
        /**
         * @brief Constructor. 
         * 
         * @param data The raw data array.
         * @param offset The offset in the raw data array. 
         * @param step The number of elements to skip between each index in the raw data array. 
         * @param length The total number of elements this Slice can access. 
         */
        MutableColumn(std::vector<T>& data, unsigned int N, unsigned int M, unsigned int col) : MutableSlice<T>(data, col, M, N) {}
        using MutableSlice<T>::MutableSlice;

        // We have to explicitly define this to avoid ambiguity
        MutableColumn& operator=(const MutableColumn& rhs) {MutableSlice<T>::operator=(rhs); return *this;}
        using MutableSlice<T>::operator=;

};