// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/Matrix.h>
#include <utility/Limit.h>
#include <io/ExistingFile.h>

namespace ausaxs {
    /**
     * @brief A representation of a dataset. The set consists of fixed number of named columns, with a variable number of rows. 
     */
    class Dataset {
        public: 
            Dataset();
            Dataset(const Dataset& d);
            Dataset(Dataset&& d);
            Dataset& operator=(const Dataset& other);
            Dataset& operator=(Dataset&& other);
            virtual ~Dataset();

            Dataset(Matrix<double>&& m);

            /**
             * @brief Create a new dataset with the given columns.
             */
            Dataset(const std::vector<std::vector<double>>& cols);

            /**
             * @brief Create a new dataset with the specified dimensions. 
             */
            Dataset(unsigned int rows, unsigned int cols);

            /**
             * @brief Create a new dataset from a data file.
             */
            Dataset(const io::ExistingFile& path);

            /**
             * @brief Get a column based on its index.
             */
            [[nodiscard]] MutableColumn<double> col(unsigned int index);

            /**
             * @brief Get a column based on its index.
             */
            [[nodiscard]] const ConstColumn<double> col(unsigned int index) const;

            /**
             * @brief Get a row based on its index.
             */
            [[nodiscard]] MutableRow<double> row(unsigned int index);

            /**
             * @brief Get a row based on its index.
             */
            [[nodiscard]] const ConstRow<double> row(unsigned int index) const;

            /**
             * @brief Get the number of points in the dataset.
             */
            [[nodiscard]] unsigned int size() const noexcept;

            /**
             * @brief Get the number of rows in the dataset.
             */
            [[nodiscard]] unsigned int size_rows() const noexcept;

            /**
             * @brief Get the number of columns in the dataset.
             */
            [[nodiscard]] unsigned int size_cols() const noexcept;

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
            Dataset select_columns(const std::vector<unsigned int>& cols) const;

            /**
             * @brief Interpolate @a num points between each pair of points in the dataset.
             */
            [[nodiscard]] Dataset interpolate(unsigned int num) const;

            /**
             * @brief Interpolate points to match the given x-values.
             */
            [[nodiscard]] Dataset interpolate(const std::vector<double>& newx) const;

            /**
             * @brief Get the interpolated value of column @a col for the given x-value.
             *        This is not suitable for looping. Use instead the interpolate(const std::vector<double>&) method.
             */
            [[nodiscard]] double interpolate_x(double x, unsigned int col) const;

            /**
             * @brief Get the weighted rolling average of this dataset. 
             *        The weight is defined as 1/(2)^i, where i is the index distance from the middle.
             * 
             * @param window The window size. 
             * 
             * @return A new (x, y) dataset with the rolling average. 
             */
            [[nodiscard]] Dataset rolling_average(unsigned int window) const;

            /**
             * @brief Get the entry with the smallest value in column @a col.
             */
            std::vector<double> find_minimum(unsigned int col) const;

            /**
             * @brief Find the indices of minima in the dataset.
             * 
             * @param min_spacing The minimum spacing between minima.
             * @param prominence The minimum prominence of a minima as a percentage of the largest prominence. Higher values will result in fewer minima.
             */
            std::vector<unsigned int> find_minima(unsigned int min_spacing = 0, double prominence = 0) const;

            /**
             * @brief Find the indices of minima in the dataset.
             * 
             * @param min_spacing The minimum spacing between minima.
             * @param prominence The minimum prominence of a minima as a percentage of the largest prominence. Higher values will result in fewer minima. 
             */
            std::vector<unsigned int> find_maxima(unsigned int min_spacing = 0, double prominence = 0) const;

            /**
             * @brief Append another dataset with the same number of rows to this one.
             *        Note that you cannot append a datasaet to itself.
             */
            void append(const Dataset& other);

            /**
             * @brief Impose limits on the data. All points with an x-value outside this range will be removed. 
             *        This assumes that the x-values are sorted. 
             *        Complexity: O(n)
             */
            void limit_x(const Limit& limits);

            /**
             * @brief Impose limits on the data. All points with an x-value outside this range will be removed. 
             *        This assumes that the x-values are sorted. 
             *        Complexity: O(n)
             */
            void limit_x(double min, double max);

            /**
             * @brief Impose limits on the data. All points with an y-value outside this range will be removed. 
             *        Complexity: O(n)
             */
            void limit_y(const Limit& limits);

            /**
             * @brief Impose limits on the data. All points with an y-value outside this range will be removed. 
             *        Complexity: O(n)
             */
            void limit_y(double min, double max);

            /**
             * @brief Sort this dataset by the x-values. 
             */
            void sort_x();

            /**
             * @brief Get the ith value in the dataset.
             */
            [[nodiscard]] double index(unsigned int i, unsigned int j) const;
            [[nodiscard]] double& index(unsigned int i, unsigned int j); //< @copydoc index(unsigned int, unsigned int) const

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
            [[nodiscard]] const ConstColumn<double> x() const {return col(0);}

            // Get the first column.
            [[nodiscard]] MutableColumn<double> x() {return col(0);}

            // Get the ith value in the first column.
            [[nodiscard]] const double& x(unsigned int i) const {return data.index(i, 0);}

            // Get the ith value in the first column.
            [[nodiscard]] double& x(unsigned int i) {return data.index(i, 0);}

            // Get the ith value in the second column.
            [[nodiscard]] const ConstColumn<double> y() const {return col(1);}

            // Get the ith value in the second column.
            [[nodiscard]] MutableColumn<double> y() {return col(1);}

            // Get the ith value in the second column.
            [[nodiscard]] const double& y(unsigned int i) const {return data.index(i, 1);}

            // Get the ith value in the second column.
            [[nodiscard]] double& y(unsigned int i) {return data.index(i, 1);}

            /**
             * @brief Define default column names.
             */
            void set_default_names();

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