#pragma once

#include <plots/PlotOptions.h>
#include <dataset/PointSet.h>
#include <dataset/Dataset.h>

namespace hist {class Histogram;}

/**
 * @brief A simple dataset is a collection of points of the form x | y | yerr. 
 */
class SimpleDataset : public Dataset, public plots::Plottable {
    protected: 
        /**
         * @brief Construct a dataset with N rows and M columns. 
         *        This is protected because it should only be used by derived classes for supporting more columns.
         */
        SimpleDataset(unsigned int N, unsigned int M);

    public: 
        /**
         * @brief Construct a new empty dataset.
         */
        SimpleDataset() noexcept;

        SimpleDataset(const hist::Histogram& h);

        SimpleDataset(hist::Histogram&& h);

        /**
         * @brief Copy constructor.
         */
        SimpleDataset(const SimpleDataset& d);

        /**
         * @brief Move constructor.
         */
        SimpleDataset(SimpleDataset&& d);

        /**
         * @brief Convert a Dataset to a SimpleDataset.
         */
        SimpleDataset(const Dataset& d);

        /**
         * @brief Construct a new empty dataset with the given number of rows. 
         */
        SimpleDataset(unsigned int rows) noexcept;

        /**
         * @brief Construct a new dataset based on the given vectors. 
         */
        SimpleDataset(std::vector<double> x, std::vector<double> y, std::vector<double> yerr);

        /**
         * @brief Construct a new dataset based on the given vectors. The errors will be initialized to 0. 
         */
        SimpleDataset(std::vector<double> x, std::vector<double> y);

        /**
         * @brief Construct a new dataset based on the given vectors. The errors will be initialized to 0. 
         */
        SimpleDataset(std::vector<double> x, std::vector<double> y, std::string xlabel, std::string ylabel);

        /**
         * @brief Construct a new dataset from an input file.
         */
        SimpleDataset(const io::ExistingFile& path);

        /**
         * @brief Destructor.
         */
        ~SimpleDataset() override = default;

        // Get the third column.
        [[nodiscard]] const ConstColumn<double> yerr() const {return col(2);}

        // Get the third column.
        [[nodiscard]] MutableColumn<double> yerr() {return col(2);}

        // Get the ith value in the third column.
        [[nodiscard]] const double& yerr(unsigned int i) const {return index(i, 2);}

        // Get the ith value in the third column.
        [[nodiscard]] double& yerr(unsigned int i) {return index(i, 2);}

        /**
         * @brief Load a dataset from the specified file. 
         */
        virtual void load(const io::ExistingFile& path) override;

        /**
         * @brief Reduce the number of rows to the specified amount by uniformly removing points in x-space.
         * 
         * @param target The target number of points.
         * @param log If true, the points will be removed uniformly on a logarithmic scale.
         */
        void reduce(unsigned int target, bool log = false);

        /**
         * @brief Assign a Matrix to this dataset.
         */
        void operator=(Matrix<double>&& other);
        
        bool operator==(const SimpleDataset& other) const;

        /**
         * @brief Get the spanned x-range. 
         */
        [[nodiscard]] Limit span_x() const noexcept;

        /**
         * @brief Get the spanned y-range. 
         */
        [[nodiscard]] Limit span_y() const noexcept;

        /**
         * @brief Get the spanned x-range.
         */
        [[nodiscard]] Limit get_xlimits() const noexcept;

        /**
         * @brief Get the spanned y-range.
         */
        [[nodiscard]] Limit get_ylimits() const noexcept;

        /**
         * @brief Get the positive spanned y-range.
         *        This can be useful for setting log ranges. 
         */
        [[nodiscard]] Limit span_y_positive() const noexcept;

        using Dataset::push_back;

        /**
         * @brief Add a new point at the end of the dataset.
         */
        void push_back(double x, double y, double yerr = 0);

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
        Point2D get_point(unsigned int index) const;

        /**
         * @brief Get the point with the smallest y-value.
         */
        Point2D find_minimum() const;

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
         * x-values will span from 0 to size-1.
         * y-values will be generated in the range [min, max]. 
         * yerr-values will be generated in the range 0.1[min, max].
         * 
         * @param size Size of the dataset.
         * @param min Minimum generated value.
         * @param max Maxium generated value. 
         */
        static SimpleDataset generate_random_data(unsigned int size, double min, double max);

        /**
         * @brief Generate a randomized dataset.
         * 
         * x-values will span from 0 to size-1.
         * y-values will be generated in the range [-value, value]. 
         * yerr-values will be generated in the range 0.1[-value value].
         * 
         * @param size Size of the dataset.
         * @param val Maximum and minimum bound on the generated values. 
         */
        static SimpleDataset generate_random_data(unsigned int size, double val);

        /**
         * @brief Get the mean of the y values.
         */
        [[nodiscard]] double mean() const;

        /**
         * @brief Get the weighted mean of the y values.
         */
        [[nodiscard]] double weighted_mean() const;

        /**
         * @brief Get the standard deviation of the y values.
         */
        [[nodiscard]] double std() const;

        /**
         * @brief Get the weighted standard deviation of the y values.
         */
        [[nodiscard]] double weighted_mean_error() const;

        /**
         * @brief Removes consecutive duplicate y-values.
         */
        void remove_consecutive_duplicates();

        using Dataset::operator=;
        SimpleDataset& operator=(const SimpleDataset& other);
        SimpleDataset& operator=(SimpleDataset&& other) noexcept;
};