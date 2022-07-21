#include <math/Matrix.h>
#include <plots/PlotOptions.h>
#include <utility/Exceptions.h>
#include <utility/PointSet.h>

class IDataset : public Matrix<double>, protected plots::PlotOptionWrapper {
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
        [[nodiscard]] size_t size() const noexcept {return N;}

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
         * @brief Get the spanned y-range. 
         */
        [[nodiscard]] Limit span_y() const noexcept;

        /**
         * @brief Get the positive spanned y-range.
         *        This can be useful for setting log ranges. 
         */
        [[nodiscard]] Limit span_y_positive() const noexcept;
};

/**
 * @brief A dataset is a collection of x and y coordinates, along with their associated errors xerr and yerr. 
 */
class Dataset : IDataset {
    public: 
        /**
         * @brief Default constructor.
         */
        Dataset() noexcept : IDataset() {
            initialize();
        }

        /**
         * @brief Construct a new empty dataset with the given number of rows.
         */
        Dataset(unsigned int rows) noexcept : IDataset(rows, 4) {
            initialize();
        }

        /**
         * @brief Construct a new dataset with x, y, xerr, and yerr values.
         */
        Dataset(std::vector<double> x, std::vector<double> y, std::vector<double> xerr, std::vector<double> yerr) noexcept : Dataset(x.size()) {
            for (unsigned int i = 0; i < x.size(); i++) {
                row(i) = {x[i], y[i], xerr[i], yerr[i]};
            }
            initialize();
        }

        /**
         * @brief Destructor.
         */
        ~Dataset() override = default;

        const ConstColumn<double> xerr() const {return ConstColumn<double>(data, N, M, 3);}
        Column<double> xerr() {return Column<double>(data, N, M, 3);}
        const double& xerr(unsigned int i) const {return index(i, 3);}
        double& xerr(unsigned int i) {return index(i, 3);}

    private:
        /**
         * @brief Initialize a dataset. 
         */
        void initialize() noexcept {
            options.xlabel = "x";
            options.ylabel = "y";
        }
};

/**
 * @brief A SAXS dataset is a collection of x and y coordinates, along with the associated error yerr. 
 */
class SAXSDataset : IDataset {
    public: 
        /**
         * @brief Construct a new empty dataset.
         */
        SAXSDataset() noexcept : IDataset() {
            initialize();
        }

        /**
         * @brief Construct a new empty dataset with the given number of rows. 
         */
        SAXSDataset(unsigned int rows) noexcept : IDataset(rows, 4) {
            initialize();
        }

        /**
         * @brief Construct a new dataset with x, y, and yerr values. 
         */
        SAXSDataset(std::vector<double> x, std::vector<double> y, std::vector<double> yerr) : SAXSDataset(x.size()) {
            for (unsigned int i = 0; i < x.size(); i++) {
                row(i) = {x[i], y[i], yerr[i]};
            }
            initialize();
        }

        /**
         * @brief Destructor.
         */
        ~SAXSDataset() override = default;

        using IDataset::push_back;

        /**
         * @brief Add a new point at the end of the dataset.
         */
        void push_back(double x, double y, double yerr) {
            extend(1);
            row(N) = {x, y, yerr};
        }

        /**
         * @brief Set the normalization of the y-values. The first y-value will be fixed to this. 
         */
        void normalize(double y0);

        /**
         * @brief Scale all errors by some common factor. 
         */
        void scale_errors(double factor);

        /**
         * @brief Scale the y-values (and their associated errors) by some common factor.
         */
        void scale_y(double factor);

        /**
         * @brief Simulate Gaussian noise on the y-values based on the errors. 
         */
        void simulate_noise();

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
         * @brief Check if this dataset has errors on the x-values.
         */
        [[nodiscard]] bool has_xerr() const noexcept;

        /**
         * @brief Check if this dataset has errors on the y-values.
         */
        [[nodiscard]] bool has_yerr() const noexcept;

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
        static SAXSDataset generate_random_data(unsigned int size, double min = 0, double max = 1);

    private: 
        /**
         * @brief Initialize a dataset.
         */
        void initialize() noexcept {
            options.xlabel = "q";
            options.ylabel = "I";
        }
};