#pragma once

#include <vector>
#include <string>
#include <memory>
#include <any>
#include <fstream>
#include <random>
#include <iterator>

#include <TGraph.h>

#include <utility/Axis.h>
#include <plots/PlotOptions.h>
#include <utility/PointSet.h>
#include <utility/View.h>

template<class T>
class IDataset {
    public: 
        /**
         * @brief Destructor. 
         */
        virtual ~IDataset() = default;

        /**
         * @brief Reduce the number of data points to the specified amount. 
         * 
         * @return The modified dataset. 
         */
        void reduce(unsigned int target, bool log = false);

        /**
         * @brief Set the plot options for this dataset.
         */
        void set_plot_options(const plots::PlotOptions& options) {plot_options = options;}

        /**
         * @brief Add to the plot options for this dataset.
         *        Any number of options can be given. 
         */
        void add_plot_options(const std::map<std::string, std::any>& options) {plot_options.set(options);}

        /**
         * @brief Add to the plot options for this dataset.
         *        This overload forces the specified style. 
         */
        void add_plot_options(std::string style, std::map<std::string, std::any> options) {plot_options.set(style, options);}

        /**
         * @brief Add to the plot options for this dataset.
         *        This overload forces the specified color. 
         */
        void add_plot_options(int color, std::map<std::string, std::any> options) {plot_options.set(color, options);}

        /**
         * @brief Set the plot color.
         */
        void set_plot_color(int color) {plot_options.color = color;}

        /**
         * @brief Get the plot options for this dataset.
         */
        virtual plots::PlotOptions get_plot_options() const {return plot_options;}

        // Get the dimension of the dataset.
        unsigned int dim() const noexcept {return T::dim();}

        // Read-only iterator to the beginning of the data.
        typename std::vector<T>::const_iterator begin() const noexcept {return data.begin();}

        // Read-only iterator to the end of the data.
        typename std::vector<T>::const_iterator end() const noexcept {return data.end();}

        // Read/write iterator to the beginning of the data.
        typename std::vector<T>::iterator begin() noexcept {return data.begin();}

        // Read/write iterator to the end of the data.
        typename std::vector<T>::iterator end() noexcept {return data.end();}

        // Push a new point to the end of the data.
        void push_back(const T& point) {data.push_back(point);}

        // Read-only indexing
        const T& operator[](unsigned int i) const {return data[i];}
        
        // Read/write indexing
        T& operator[](unsigned int i) {return data[i];}

        // Get the number of points in the dataset.
        unsigned int size() const noexcept {return data.size();}

        std::vector<T> data;

    protected: 
        plots::PlotOptions plot_options;
};

namespace detail {
    template<class T>
    class XOps : public virtual IDataset<T> {
        static_assert(std::is_base_of<detail::IPoint, T>::value, "T must be a subclass of IPoint");
        static_assert(T::dim() >= 1, "X operations are only supported with 1D data or higher.");

        public:
            /**
             * @brief Get a view of the x column.
             */
            XView<T> x() {return XView<T>(this->data);}

            /**
             * @brief Get a view of the x column.
             */
            const XConstView<T> x() const {return XConstView<T>(this->data);}

            /**
             * @brief Get a specific x value. 
             */
            double& x(unsigned int index) {return this->data[index].x;}

            /**
             * @brief Get a specific x value. 
             */
            const double& x(unsigned int index) const {return this->data[index].x;}

            /**
             * @brief Impose a limit on the x values.
             */
            void limit_x(const Limit& limit) {
                std::vector<T> new_data;
                for (auto& point : this->data) {
                    if (point.x < limit.min) {continue;} 
                    else if (point.x > limit.max) {continue;}
                    new_data.push_back(point);
                }
                this->data = std::move(new_data);
            }

            /**
             * @brief Get the bounds on the x values.
             */
            Limit get_limit_x() const {
                Limit limit;
                for (const auto& point : this->data) {
                    if (point.x < limit.min) {limit.min = point.x;}
                    else if (point.x > limit.max) {limit.max = point.x;}
                }
                return limit;
            }

        private: 
            class XIterator {
                using value_type = T;
                using difference_type = std::ptrdiff_t;
                using pointer = T*;
                using reference = T&;
                using iterator_category = std::forward_iterator_tag;

                public: 
                    XIterator(unsigned int index) : index(index) {}

                    reference operator*() override {return this->data[this->index].x;}
                    pointer operator->() override {return &this->data[this->index].x;}

                    XIterator& operator++() {++index; return *this;}
                    XIterator operator++(int) {XIterator tmp(*this); ++(*this); return tmp;}

                    XIterator& operator--() {--index; return *this;}
                    XIterator operator--(int) {XIterator tmp(*this); --(*this); return tmp;}

                    friend bool operator==(const XIterator& lhs, const XIterator& rhs) {return lhs.index == rhs.index;}
                    friend bool operator!=(const XIterator& lhs, const XIterator& rhs) {return !(lhs.index == rhs.index);}
                
                private:
                    unsigned int index;
            };

            // Read/write iterator to the beginning of the data.
            XIterator begin() noexcept {return XIterator();}

            // Read/write iterator to the end of the data.
            XIterator end() noexcept {return XIterator();}
    };

    template<class T>
    class XErrOps : public virtual IDataset<T> {
        static_assert(std::is_base_of<detail::IPointError, T>::value, "T must be a subclass of IPoint");
        static_assert(T::dim() >= 1, "Xerr operations are only supported with 1D data or higher.");

        public:
            /**
             * @brief Get a view of the xerr column.
             */
            XErrorView<T> xerr() {return XErrorView<T>(this->data);}

            /**
             * @brief Get a view of the xerr column.
             */
            const XConstErrorView<T> xerr() const {return XConstErrorView<T>(this->data);}

            /**
             * @brief Get a specific xerr value. 
             */
            double& xerr(unsigned int index) {return this->data[index].xerr;}

            /**
             * @brief Get a specific xerr value. 
             */
            const double& xerr(unsigned int index) const {return this->data[index].xerr;}
    };

    template<class T>
    class YOps : public virtual IDataset<T> {
        static_assert(std::is_base_of<detail::IPoint, T>::value, "T must be a subclass of IPoint");
        static_assert(T::dim() >= 2, "Y operations are only supported with 2D data or higher.");

        public:
            /**
             * @brief Get a view of the y column.
             */
            YView<T> y() {return YView<T>(this->data);}

            /**
             * @brief Get a view of the y column.
             */
            const YConstView<T> y() const {return YConstView<T>(this->data);}

            /**
             * @brief Get a specific y value. 
             */
            double& y(unsigned int index) {return this->data[index].y;}

            /**
             * @brief Get a specific y value. 
             */
            const double& y(unsigned int index) const {return this->data[index].y;}

            /**
             * @brief Impose a limit on the x values.
             */
            void limit_y(const Limit& limit) {
                std::vector<T> new_data;
                for (auto& point : this->data) {
                    if (point.y < limit.min) {continue;} 
                    else if (point.y > limit.max) {continue;}
                    new_data.push_back(point);
                }
                this->data = std::move(new_data);
            }

            /**
             * @brief Get the bounds on the x values.
             */
            Limit get_limit_y() const {
                Limit limit;
                for (const auto& point : this->data) {
                    if (point.y < limit.min) {limit.min = point.y;}
                    else if (point.y > limit.max) {limit.max = point.y;}
                }
                return limit;
            }
    };

    template<class T>
    class YErrOps : public virtual IDataset<T> {
        static_assert(std::is_base_of<detail::IPointError, T>::value, "T must be a subclass of IPoint");
        static_assert(T::dim() >= 2, "Yerr operations are only supported with 2D data or higher.");

        public:
            /**
             * @brief Get a view of the yerr column.
             */
            YErrorView<T> yerr() {return YErrorView<T>(this->data);}

            /**
             * @brief Get a view of the yerr column.
             */
            const YConstErrorView<T> yerr() const {return YConstErrorView<T>(this->data);}

            /**
             * @brief Get a specific yerr value. 
             */
            double& yerr(unsigned int index) {return this->data[index].yerr;}

            /**
             * @brief Get a specific yerr value. 
             */
            const double& yerr(unsigned int index) const {return this->data[index].yerr;}
    };
}

/**
 * @brief A collection of points with x- and y-coordinates.
 */
class PointSet : virtual public detail::XOps<Point2D>, virtual public detail::YOps<Point2D> {
    public:
        PointSet() noexcept = default;

        virtual ~PointSet() override = default;
};

/**
 * @brief A collection of data points with x- and y-coordinates and associated errors.
 */
class Dataset : virtual public detail::XOps<ErrorPoint2D>, virtual public detail::YOps<ErrorPoint2D>, virtual public detail::XErrOps<ErrorPoint2D>, virtual public detail::YErrOps<ErrorPoint2D> {
    public:
        /**
         * @brief Constructor. 
         */
        Dataset() = default;

        /**
         * @brief Construct a dataset from a vector of points and associated errors.
         */
        Dataset(std::vector<double> x, std::vector<double> y, std::vector<double> xerr, std::vector<double> yerr) noexcept {
            data.resize(x.size());
            for (unsigned int i = 0; i < x.size(); i++) {
                data[i].x = x[i];
                data[i].y = y[i];
                data[i].xerr = xerr[i];
                data[i].yerr = yerr[i];
            }
        }

        Dataset(std::vector<double> x, std::vector<double> y, std::string xlabel, std::string ylabel) : Dataset(x, y, std::vector<double>(x.size(), 0.0), std::vector<double>(y.size(), 0.0)) {
            plot_options.xlabel = xlabel;
            plot_options.ylabel = ylabel;
        }

        /**
         * @brief Construct a dataset with x, y, and yerr columns. The xerr column is initialized to 0.
         */
        Dataset(std::vector<double> x, std::vector<double> y, std::vector<double> yerr) noexcept : Dataset(x, y, std::vector<double>(x.size(), 0), yerr) {}

        /**
         * @brief Construct a dataset with x, y columns. The error columns are initialized to 0.
         */
        Dataset(std::vector<double> x, std::vector<double> y) noexcept : Dataset(x, y, std::vector<double>(x.size(), 0), std::vector<double>(x.size(), 0)) {}

        /**
         * @brief Destructor. 
         */
        virtual ~Dataset() override = default;

        /**
         * @brief Check if this dataset contains y errors.
         */
        bool has_yerr() const {
            if (data.size() == 0) {return false;}
            return data[0].yerr != 0;
        }

        /**
         * @brief Check if this dataset contains x errors.
         */
        bool has_xerr() const {
            if (data.size() == 0) {return false;}
            return data[0].xerr != 0;
        }

        /**
         * @brief Add a point at the end of the dataset.
         */
        void push_back(double x, double y) {
            data.push_back(ErrorPoint2D(x, y, 0, 0));
        }

        /**
         * @brief Add a point at the end of the dataset.
         */
        void push_back(double x, double y, double xerr, double yerr) {
            data.push_back(ErrorPoint2D(x, y, xerr, yerr));
        }

        /**
         * @brief Simulate Gaussian noise on the y-values based on the errors. 
         */
        void simulate_noise();

        /**
         * @brief Generate errors for the y-values mimicking what one would find experimentally. 
         */
        void simulate_errors();

        /**
         * @brief Set the resolution of this dataset. 
         */
        void set_resolution(unsigned int resolution);

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
        static Dataset generate_random_data(unsigned int size, double min = 0, double max = 1);

        /**
         * @brief Reduce the number of data points to the specified amount. 
         * 
         * @return The modified dataset. 
         */
        void reduce(unsigned int target, bool log = false);

        /**
         * @brief Get the point with the smallest y-value.
         */
        Point2D find_minimum() const noexcept;

        /**
         * @brief Get a set of data based on its label. 
         * This also serves as a consistency check. 
         * 
         * @param label The name of the data. 
         */
        std::vector<double>& get(const std::string label);

        /**
         * @brief Check if the data is logarithmic. 
         *        This may be wrong if the x-data is very noisy.
         */
        bool is_logarithmic() const noexcept;

        /**
         * @brief Scale all errors by some common factor. 
         */
        void scale_errors(double factor);

        /**
         * @brief Scale the y-values (and their associated errors) by some common factor.
         */
        void scale_y(double factor);

        /**
         * @brief Set the normalization of the y-values. The first y-value will be fixed to this. 
         */
        void normalize(double y0);

        /**
         * @brief Plot this dataset.
         */
        std::unique_ptr<TGraph> plot() const;

        /**
         * @brief Get the spanned y-range. 
         */
        [[nodiscard]] Limit span() const noexcept;

        /**
         * @brief Get the positive spanned y-range.
         *        This can be useful for setting log ranges. 
         */
        [[nodiscard]] Limit span_positive() const noexcept;

        /**
         * @brief Write this dataset to the specified file. 
         */
        void save(std::string path) const;

        /**
         * @brief Get the plot options for this dataset.
         */
        virtual plots::PlotOptions get_plot_options() const override {
            plots::PlotOptions options = plot_options;
            if (options.xlabel.empty()) {
                options.xlabel = "x";
            }
            if (options.ylabel.empty()) {
                options.ylabel = "y";
            }
            return options;
        }
};

// Define a specialized SAXSDataset class for additional type safety. This makes it very clear what this dataset is supposed to contain.
class SAXSDataset : public Dataset {
    using Dataset::Dataset;

    public:
        /**
         * @brief Construct a dataset from an input file.
         */
        SAXSDataset(std::string file) {
            load(file);
        }

    private: 
        /**
         * @brief Load a dataset from a file.
         * 
         * @param file Path to the input file. 
         */
        void load(const std::string file);

        /**
         * @brief Get the plot options for this dataset.
         */
        virtual plots::PlotOptions get_plot_options() const override {
            plots::PlotOptions options = plot_options;
            if (options.xlabel.empty()) {
                options.xlabel = "q";
            }
            if (options.ylabel.empty()) {
                options.ylabel = "I";
            }
            return options;
        }
};