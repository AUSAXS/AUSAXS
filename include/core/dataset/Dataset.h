// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <io/ExistingFile.h>
#include <math/Matrix.h>
#include <utility/Limit.h>

namespace ausaxs {
    /**
     * @brief A representation of a dataset. The set consists of fixed number of named columns, with a variable number of rows. 
     */
    class Dataset {
        public: 
            Dataset();
            Dataset(const Dataset& d);
            Dataset(Dataset&& d) noexcept ;
            Dataset& operator=(const Dataset& other);
            Dataset& operator=(Dataset&& other) noexcept ;
            virtual ~Dataset();

            template<typename ...Args> requires std::constructible_from<Matrix<double>, Args...>
            Dataset(Args&&... args) : data(std::forward<Args>(args)...) {}
            Dataset(std::initializer_list<std::vector<double>> lists) : data(lists) {}

            /**
             * @brief Create a new dataset from a data file.
             */
            Dataset(const io::ExistingFile& path);

            /**
             * @brief Get a column based on its index.
             */
            [[nodiscard]] MutableColumn<double> col(int index);

            /**
             * @brief Get a column based on its index.
             */
            [[nodiscard]] ConstColumn<double> col(int index) const;

            /**
             * @brief Get a row based on its index.
             */
            [[nodiscard]] MutableRow<double> row(int index);

            /**
             * @brief Get a row based on its index.
             */
            [[nodiscard]] ConstRow<double> row(int index) const;

            /**
             * @brief Get the number of points in the dataset.
             */
            [[nodiscard]] int size() const noexcept;

            /**
             * @brief Get the number of rows in the dataset.
             */
            [[nodiscard]] int size_rows() const noexcept;

            /**
             * @brief Get the number of columns in the dataset.
             */
            [[nodiscard]] int size_cols() const noexcept;

            /**
             * @brief Check if the dataset is empty.
             */
            [[nodiscard]] bool empty() const noexcept;

            /**
             * @brief Write this dataset to the specified file. 
             * 
             * @param path The path to the save location.
             * @param header The header for the file. 
             */
            void save(const io::File& path, const std::string& header = "") const;

            /**
            * @brief Create a new dataset with the specified columns.
            */
            Dataset select_columns(const std::vector<int>& cols) const;

            /**
             * @brief Interpolate @a num points between each pair of points in the dataset.
             */
            [[nodiscard]] Dataset interpolate(int num) const;

            /**
             * @brief Interpolate points to match the given x-values.
             */
            [[nodiscard]] Dataset interpolate(const std::vector<double>& newx) const;

            /**
             * @brief Get the interpolated value of column @a col for the given x-value.
             *        This is not suitable for looping. Use instead the interpolate(const std::vector<double>&) method.
             */
            [[nodiscard]] double interpolate_x(double x, int col) const;

            /**
             * @brief Get the weighted rolling average of this dataset. 
             *        The weight is defined as 1/(2)^i, where i is the index distance from the middle.
             * 
             * @param window The window size. 
             * 
             * @return A new (x, y) dataset with the rolling average. 
             */
            [[nodiscard]] Dataset rolling_average(int window) const;

            /**
             * @brief Get the entry with the smallest value in column @a col.
             */
            std::vector<double> find_minimum(int col) const;

            /**
             * @brief Get the range spanned by the values in column @a col.
             */
            [[nodiscard]] Limit span(int col) const noexcept;

            /**
             * @brief Get the mean of the values in column @a col.
             */
            [[nodiscard]] double mean(int col) const;

            /**
             * @brief Get the standard deviation of the values in column @a col.
             */
            [[nodiscard]] double std(int col) const;

            /**
             * @brief Append another dataset with the same number of rows to this one.
             *        Note that you cannot append a datasaet to itself.
             */
            void append(const Dataset& other);

            /**
             * @brief Impose limits on the data. All rows with a value in column @a col outside this range will be removed. 
             *        Complexity: O(n)
             */
            void limit(int col, const Limit& limits);

            /**
             * @brief Impose limits on the data. All rows with a value in column @a col outside this range will be removed. 
             *        Complexity: O(n)
             */
            void limit(int col, double min, double max);

            /**
             * @brief Sort the rows of this dataset by the values in column @a col. 
             */
            void sort(int col);

            /**
             * @brief Get the ith value in the dataset.
             */
            [[nodiscard]] double index(int i, int j) const;
            [[nodiscard]] double& index(int i, int j); //< @copydoc index(int, int) const

            /**
             * @brief Add a new row to the dataset.
             */
            void push_back(const std::vector<double>& row);

            /**
             * @brief Get the string representation of this object.
             */
            [[nodiscard]] std::string to_string() const;

            [[nodiscard]] bool operator==(const Dataset& other) const;

        //#####################//
        //### Alias methods ###//
        //#####################//

            // Get the first column.
            [[nodiscard]] ConstColumn<double> x() const {return col(0);}

            // Get the first column.
            [[nodiscard]] MutableColumn<double> x() {return col(0);}

            // Get the ith value in the first column.
            [[nodiscard]] const double& x(int i) const {return data.index(i, 0);}

            // Get the ith value in the first column.
            [[nodiscard]] double& x(int i) {return data.index(i, 0);}

            // Get the ith value in the second column.
            [[nodiscard]] ConstColumn<double> y() const {return col(1);}

            // Get the ith value in the second column.
            [[nodiscard]] MutableColumn<double> y() {return col(1);}

            // Get the ith value in the second column.
            [[nodiscard]] const double& y(int i) const {return data.index(i, 1);}

            // Get the ith value in the second column.
            [[nodiscard]] double& y(int i) {return data.index(i, 1);}

            Matrix<double> data;
 
        protected:
            /**
             * @brief Load a dataset from the specified file. 
             */
            virtual void load(const io::ExistingFile& path);

            /**
             * @brief Assign a matrix to this Dataset.
             *        The compatibility of the matrix is checked.
             */
            void assign_matrix(Matrix<double>&& m);

            /**
             * @brief Forcibly assign a matrix to this Dataset.
             *        The compatibility of the matrix is not checked. 
             */
            void force_assign_matrix(Matrix<double>&& m);
    };
}