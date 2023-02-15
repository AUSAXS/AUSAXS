#pragma once

#include <math/Matrix.h>
#include <utility/Exceptions.h>
#include <dataset/PointSet.h>

/**
 * @brief A representation of a dataset. The set consists of fixed number of named columns, with a variable number of rows. 
 */
class Dataset : public Matrix<double> {
    public: 
        /**
         * @brief Default constructor. 
         */
        Dataset() {}

        /**
         * @brief Matrix constructor.
         */
        Dataset(Matrix&& m) : Matrix(std::move(m)) {}

        /**
         * @brief Create a new dataset with the given columns.
         */
        Dataset(std::vector<std::string> col_names) : Matrix(0, col_names.size()), names(col_names) {}

        /**
         * @brief Create a new dataset with the given columns.
         */
        Dataset(std::vector<std::vector<double>> cols, std::vector<std::string> col_names) : Matrix(cols), names(col_names) {}

        /**
         * @brief Create a new dataset with the given columns.
         */
        Dataset(std::vector<std::vector<double>> cols);

        /**
         * @brief Create a new dataset with the specified dimensions. 
         */
        Dataset(unsigned int rows, unsigned int cols);

        /**
         * @brief Destructor.
         */
        virtual ~Dataset() = default;

        /**
         * @brief Get a column based on its name. 
         */
        [[nodiscard]] Column<double> col(std::string column);

        /**
         * @brief Get a column based on its name. 
         */
        [[nodiscard]] const ConstColumn<double> col(std::string column) const;

        /**
         * @brief Get a column based on its index.
         */
        [[nodiscard]] Column<double> col(unsigned int index);

        /**
         * @brief Get a column based on its index.
         */
        [[nodiscard]] const ConstColumn<double> col(unsigned int index) const;

        /**
         * @brief Get a row based on its index.
         */
        [[nodiscard]] Row<double> row(unsigned int index);

        /**
         * @brief Get a row based on its index.
         */
        [[nodiscard]] const ConstRow<double> row(unsigned int index) const;

        /**
         * @brief Get the number of points in the dataset.
         */
        [[nodiscard]] size_t size() const noexcept;

        [[nodiscard]] bool empty() const noexcept;

        /**
         * @brief Write this dataset to the specified file. 
         * 
         * @param path The path to the save location.
         * @param header The header for the file. 
         */
        void save(std::string path, std::string header = "") const;

        /**
         * @brief Load a dataset from the specified file. 
         */
        virtual void load(std::string path);

        /**
         * @brief Set the column names. 
         */
        void set_col_names(std::vector<std::string> names);

        /**
         * @brief Set a column name. 
         */
        void set_col_names(unsigned int i, std::string name);

        /**
         * @brief Get the column names. 
         */
        [[nodiscard]] std::vector<std::string> get_col_names();

        /**
         * @brief Get a column name. 
         */
        [[nodiscard]] std::string get_col_names(unsigned int i);

        /**
         * @brief Interpolate @a num points between each pair of points in the dataset.
         */
        void interpolate(unsigned int num);

        /**
         * @brief Get the weighted rolling average of this dataset. 
         *        The weight is defined as 1/(2)^i, where i is the index distance from the middle.
         * 
         * @param window The window size. 
         * 
         * @return A new (x, y) dataset with the rolling average. 
         */
        Dataset rolling_average(unsigned int window) const;

        /**
         * @brief Append another dataset to this one.
         */
        void append(const Dataset& other);

        std::string to_string() const;

    //*#####################*//
    //*### Alias methods ###*//
    //*#####################*//

        // Get the first column.
        [[nodiscard]] const ConstColumn<double> x() const {return col(0);}

        // Get the first column.
        [[nodiscard]] Column<double> x() {return col(0);}

        // Get the ith value in the first column.
        [[nodiscard]] const double& x(unsigned int i) const {return index(i, 0);}

        // Get the ith value in the first column.
        [[nodiscard]] double& x(unsigned int i) {return index(i, 0);}

        // Get the ith value in the second column.
        [[nodiscard]] const ConstColumn<double> y() const {return col(1);}

        // Get the ith value in the second column.
        [[nodiscard]] Column<double> y() {return col(1);}

        // Get the ith value in the second column.
        [[nodiscard]] const double& y(unsigned int i) const {return index(i, 1);}

        // Get the ith value in the second column.
        [[nodiscard]] double& y(unsigned int i) {return index(i, 1);}

        std::vector<std::string> names; // The column names
    private: 
        /**
         * @brief Define default column names.
         */
        void set_default_names();

    protected:
        /**
         * @brief Assign a matrix to this Dataset.
         */
        void assign_matrix(const Matrix<double>&& m);
};