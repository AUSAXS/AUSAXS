#pragma once

#include <math/Matrix.h>
#include <utility/Exceptions.h>
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
         * @brief Create a new dataset with the given columns.
         */
        Dataset(std::vector<std::string> col_names) : Matrix(0, col_names.size()), names(col_names) {}

        /**
         * @brief Create a new dataset with the given columns.
         */
        Dataset(std::vector<std::vector<double>> cols) : Matrix(cols) {
            for (unsigned int i = 0; i < cols.size(); i++) {
                names.push_back(std::to_string(i));
            }
        }

        /**
         * @brief Create a new dataset with the given columns.
         */
        Dataset(std::vector<std::vector<double>> cols, std::vector<std::string> col_names) : Matrix(cols), names(col_names) {}

        /**
         * @brief Create a new dataset with the specified dimensions. 
         */
        Dataset(unsigned int rows, unsigned int cols) : Matrix(rows , cols) {}

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

        /**
         * @brief Assign a Matrix to this dataset.
         */
        void operator=(const Matrix<double>&& other);

        /**
         * @brief Add a point at the end of the dataset.
         */
        void push_back(double x, double y, double xerr, double yerr) {
            data.push_back(ErrorPoint2D(x, y, xerr, yerr));
        }

        /**
         * @brief Load a dataset from the specified file. 
         */
        void load(std::string path);

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

    private: 
        /**
         * @brief Get a column name. 
         */
        [[nodiscard]] std::string get_col_names(unsigned int i);

    private: 
        std::vector<std::string> names; // The column names
};