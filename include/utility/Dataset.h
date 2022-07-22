#pragma once

#include <math/Matrix.h>
#include <plots/PlotOptions.h>
#include <utility/Exceptions.h>
#include <utility/PointSet.h>

class IDataset : public Matrix<double>, public plots::PlotOptionWrapper {
    public: 
        using Matrix::Matrix;
        virtual ~IDataset() = default;

        [[nodiscard]] const ConstColumn<double> x() const {return ConstColumn<double>(data, N, M, 0);}
        [[nodiscard]] Column<double> x() {return Column<double>(data, N, M, 0);}
        [[nodiscard]] const double& x(unsigned int i) const {return index(i, 0);}
        [[nodiscard]] double& x(unsigned int i) {return index(i, 0);}

        [[nodiscard]] const ConstColumn<double> y() const {return ConstColumn<double>(data, N, M, 1);}
        [[nodiscard]] Column<double> y() {return Column<double>(data, N, M, 1);}
        [[nodiscard]] const double& y(unsigned int i) const {return index(i, 1);}
        [[nodiscard]] double& y(unsigned int i) {return index(i, 1);}

        [[nodiscard]] const ConstColumn<double> yerr() const {return ConstColumn<double>(data, N, M, 2);}
        [[nodiscard]] Column<double> yerr() {return Column<double>(data, N, M, 2);}
        [[nodiscard]] const double& yerr(unsigned int i) const {return index(i, 2);}
        [[nodiscard]] double& yerr(unsigned int i) {return index(i, 2);}

        /**
         * @brief Get a column based on its name. This is primarily meant as a sanity check when expecting e.g. logarithmic data. 
         */
        [[nodiscard]] Column<double> get(std::string column) {
            if (column == options.xlabel) {return x();}
            else if (column == options.ylabel) {return y();}
            else {throw except::invalid_argument("Error in IDataset::get: Column name \"" + column + "\" not recognized.");}
        }

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
         *        This assumes that the y-values are unsorted. 
         *        Complexity: O(n)
         */
        void limit_y(const Limit& limits);

        /**
         * @brief Impose limits on the data. All points with an y-value outside this range will be removed. 
         *        This assumes that the y-values are unsorted. 
         *        Complexity: O(n)
         */
        void limit_y(double min, double max);

        /**
         * @brief Get the number of points in the dataset.
         */
        [[nodiscard]] size_t size() const noexcept;

        /**
         * @brief Reduce the number of data points to the specified amount. 
         */
        void reduce(unsigned int target, bool log = false);

        /**
         * @brief Assign a Matrix to this dataset.
         */
        void operator=(const Matrix<double>&& other);

        /**
         * @brief Check if the data is logarithmic. 
         *        This may be wrong if the x-data is very noisy.
         */
        bool is_logarithmic() const noexcept;

        /**
         * @brief Write this dataset to the specified file. 
         */
        void save(std::string path) const;

        /**
         * @brief Get the spanned x-range. 
         */
        [[nodiscard]] Limit span_x() const noexcept;

        /**
         * @brief Get the spanned y-range. 
         */
        [[nodiscard]] Limit span_y() const noexcept;

        /**
         * @brief Get the positive spanned y-range.
         *        This can be useful for setting log ranges. 
         */
        [[nodiscard]] Limit span_y_positive() const noexcept;

        /**
         * @brief Load a dataset from the specified file. 
         */
        void load(std::string path);
};

/**
 * @brief A simple dataset is a collection of points of the form x | y | yerr. 
 */
class SimpleDataset : public IDataset {
    protected: 
        /**
         * @brief Construct a dataset with N rows and M columns. 
         *        This is protected because it should only be used by derived classes for supporting more columns.
         */
        SimpleDataset(unsigned int N, unsigned int M) : IDataset(N, M) {}

    public: 
        /**
         * @brief Construct a new empty dataset with the given number of rows. 
         */
        SimpleDataset(unsigned int rows) noexcept : IDataset(rows, 3) {}

        /**
         * @brief Construct a new empty dataset.
         */
        SimpleDataset() noexcept : SimpleDataset(0) {}

        /**
         * @brief Construct a new dataset based on the given vectors. 
         */
        SimpleDataset(std::vector<double> x, std::vector<double> y, std::vector<double> yerr) : SimpleDataset(x.size()) {
            if (x.size() != y.size() || x.size() != yerr.size()) {
                throw except::size_error("Error in SimpleDataset::SimpleDataset: x, y, and yerr must have the same size.");
            }
            for (unsigned int i = 0; i < x.size(); i++) {
                row(i) = {x[i], y[i], yerr[i]};
            }
        }

        /**
         * @brief Construct a new dataset based on the given vectors. The errors will be initialized to 0. 
         */
        SimpleDataset(std::vector<double> x, std::vector<double> y) : SimpleDataset(x, y, std::vector<double>(x.size(), 0)) {}

        /**
         * @brief Construct a new dataset based on the given vectors. The errors will be initialized to 0. 
         */
        SimpleDataset(std::vector<double> x, std::vector<double> y, std::string xlabel, std::string ylabel) : SimpleDataset(x, y, std::vector<double>(x.size(), 0)) {
            options.xlabel = xlabel;
            options.ylabel = ylabel;
        }

        /**
         * @brief Construct a new dataset from an input file.
         */
        SimpleDataset(std::string path) {
            load(path);
        }

        /**
         * @brief Destructor.
         */
        ~SimpleDataset() override = default;

        using IDataset::push_back;

        /**
         * @brief Add a new point at the end of the dataset.
         */
        void push_back(double x, double y, double yerr) {
            extend(1);
            row(N) = {x, y, yerr};
        }

        /**
         * @brief Name the columns. 
         */
        void name_columns(std::string xlabel, std::string ylabel) {
            options.xlabel = xlabel;
            options.ylabel = ylabel;
        }

        /**
         * @brief Set the normalization of the y-values. The first y-value will be fixed to this. 
         */
        void normalize(double y0);

        /**
         * @brief Scale all errors by some common factor. 
         */
        virtual void scale_errors(double factor);

        /**
         * @brief Scale the y-values (and their associated errors) by some common factor.
         */
        void scale_y(double factor);

        /**
         * @brief Simulate Gaussian noise on the y-values based on the errors. 
         */
        void simulate_noise();

        /**
         * @brief Generate errors for the y-values mimicking what one would find experimentally. 
         */
        void simulate_errors();

        /**
         * @brief Get the point at a given index.
         */
        Point2D get_point(unsigned int index) const noexcept;

        /**
         * @brief Get the point with the smallest y-value.
         */
        Point2D find_minimum() const noexcept;

        /**
         * @brief Add a new datapoint to the end of this dataset. 
         */
        void push_back(const Point2D& point) noexcept;

        /**
         * @brief Rebin the data to a logarithmic scale. 
         *        This follows the typical rebinning algorithm used experimentally.
         */
        void rebin() noexcept;

        /**
         * @brief Generate a randomized dataset.
         * 
         * @param size Size of the dataset.
         * @param min Minimum generated value.
         * @param max Maxium generated value. 
         */
        static SimpleDataset generate_random_data(unsigned int size, double min = 0, double max = 1);
};


/**
 * @brief A dataset is a collection of points of the form x | y | xerr | yerr. 
 */
class Dataset : public SimpleDataset {
    public: 
        using SimpleDataset::SimpleDataset;

        /**
         * @brief Default constructor.
         */
        Dataset() noexcept : SimpleDataset() {}

        /**
         * @brief Construct a new empty dataset with the given number of rows.
         */
        Dataset(unsigned int rows) noexcept : SimpleDataset(rows, 4) {}

        /**
         * @brief Construct a new dataset with x, y, xerr, and yerr values.
         */
        Dataset(std::vector<double> x, std::vector<double> y, std::vector<double> xerr, std::vector<double> yerr) noexcept : Dataset(x.size()) {
            for (unsigned int i = 0; i < x.size(); i++) {
                row(i) = {x[i], y[i], xerr[i], yerr[i]};
            }
        }

        Dataset(SimpleDataset data) : Dataset(data.size()) {
            for (unsigned int i = 0; i < data.size(); i++) {
                row(i) = {data.x(i), data.y(i), data.yerr(i), 0};
            }
        }

        /**
         * @brief Destructor.
         */
        ~Dataset() override = default;

        /**
         * @brief Add a new point at the end of the dataset.
         */
        void push_back(double x, double y, double xerr, double yerr) {
            extend(1);
            row(N) = {x, y, yerr, xerr};
        }

        /**
         * @brief Add a new point at the end of the dataset.
         */
        void push_back(double x, double y) {
            push_back(x, y, 0, 0);
        }

        /**
         * @brief Add a new point at the end of the dataset.
         */
        void push_back(const Point2D& point) noexcept {
            push_back(point.x, point.y, point.xerr, point.yerr);
        }

        /**
         * @brief Scale all errors by some common factor. 
         */
        void scale_errors(double factor) override;

        [[nodiscard]] const ConstColumn<double> xerr() const {return ConstColumn<double>(data, N, M, 3);}
        [[nodiscard]] Column<double> xerr() {return Column<double>(data, N, M, 3);}
        [[nodiscard]] const double& xerr(unsigned int i) const {return index(i, 3);}
        [[nodiscard]] double& xerr(unsigned int i) {return index(i, 3);}
};

// Object conversion between Dataset and SimpleDataset is often used. This check ensures that the conversion is safe.
static_assert(sizeof(Dataset) == sizeof(SimpleDataset), "Object conversion is broken.");