#pragma once

#include <utility/PointSet.h>
#include <utility/SimpleDataset.h>

/**
 * @brief A dataset is a collection of points of the form x | y | xerr | yerr. 
 */
class Dataset2D : public SimpleDataset {
    public: 
        /**
         * @brief Default constructor.
         */
        Dataset2D() noexcept : SimpleDataset(0, 4) {}

        /**
         * @brief Construct a new empty dataset with the given number of rows.
         */
        Dataset2D(unsigned int rows) noexcept : SimpleDataset(rows, 4) {}

        /**
         * @brief Construct a new dataset with x and y values. The xerr and yerr columns will be initialized to 0.
         */
        Dataset2D(std::vector<double> x, std::vector<double> y) noexcept : Dataset2D(x.size()) {
            for (unsigned int i = 0; i < x.size(); i++) {
                row(i) = {x[i], y[i], 0, 0};
            }
        }

        /**
         * @brief Construct a new dataset based on the given vectors. The errors will be initialized to 0. 
         */
        Dataset2D(std::vector<double> x, std::vector<double> y, std::string xlabel, std::string ylabel) : Dataset2D(x, y) {
            set_col_names({xlabel, ylabel, std::string(ylabel)+"err", std::string(xlabel)+"err"});
            options.xlabel = xlabel;
            options.ylabel = ylabel;
        }

        /**
         * @brief Construct a new dataset with x, y, and yerr values. The xerr column will be initialized to 0. 
         */
        Dataset2D(std::vector<double> x, std::vector<double> y, std::vector<double> yerr) noexcept : Dataset2D(x.size()) {
            for (unsigned int i = 0; i < x.size(); i++) {
                row(i) = {x[i], y[i], yerr[i], 0};
            }
        }

        /**
         * @brief Construct a new dataset with x, y, xerr, and yerr values.
         */
        Dataset2D(std::vector<double> x, std::vector<double> y, std::vector<double> xerr, std::vector<double> yerr) noexcept : Dataset2D(x.size()) {
            for (unsigned int i = 0; i < x.size(); i++) {
                row(i) = {x[i], y[i], xerr[i], yerr[i]};
            }
        }

        Dataset2D(const SimpleDataset& data) : Dataset2D(data.size()) {
            for (unsigned int i = 0; i < data.size(); i++) {
                row(i) = {data.x(i), data.y(i), data.yerr(i), 0};
            }
        }

        /**
         * @brief Construct a new dataset from an input file.
         */
        Dataset2D(std::string path) : Dataset2D() {
            load(path);
        }

        /**
         * @brief Destructor.
         */
        ~Dataset2D() override = default;

        /**
         * @brief Add a new point at the end of the dataset.
         */
        void push_back(double x, double y, double xerr, double yerr);

        /**
         * @brief Add a new point at the end of the dataset.
         */
        void push_back(double x, double y);

        /**
         * @brief Add a new point at the end of the dataset.
         */
        void push_back(const Point2D& point) noexcept;

        /**
         * @brief Scale all errors by some common factor. 
         */
        void scale_errors(double factor) override;

        [[nodiscard]] const ConstColumn<double> xerr() const {return ConstColumn<double>(data, N, M, 3);}
        [[nodiscard]] Column<double> xerr() {return Column<double>(data, N, M, 3);}
        [[nodiscard]] const double& xerr(unsigned int i) const {return index(i, 3);}
        [[nodiscard]] double& xerr(unsigned int i) {return index(i, 3);}

    private: 
        void load(std::string path) override;
};

// Object conversion between Dataset2D and SimpleDataset is often used. This check ensures that the conversion is safe.
static_assert(sizeof(Dataset2D) == sizeof(SimpleDataset), "Object conversion is broken.");