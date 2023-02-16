#pragma once

#include <dataset/PointSet.h>
#include <dataset/SimpleDataset.h>

/**
 * @brief A dataset is a collection of points of the form x | y | xerr | yerr. 
 */
class Dataset2D : public SimpleDataset {
    public: 
        /**
         * @brief Default constructor.
         */
        Dataset2D() noexcept;

        /**
         * @brief Construct a new empty dataset with the given number of rows.
         */
        Dataset2D(unsigned int rows) noexcept;

        /**
         * @brief Construct a new dataset with x and y values. The xerr and yerr columns will be initialized to 0.
         */
        Dataset2D(std::vector<double> x, std::vector<double> y) noexcept;

        /**
         * @brief Construct a new dataset based on the given vectors. The errors will be initialized to 0. 
         */
        Dataset2D(std::vector<double> x, std::vector<double> y, std::string xlabel, std::string ylabel);

        /**
         * @brief Construct a new dataset with x, y, and yerr values. The xerr column will be initialized to 0. 
         */
        Dataset2D(std::vector<double> x, std::vector<double> y, std::vector<double> yerr) noexcept;

        /**
         * @brief Construct a new dataset with x, y, xerr, and yerr values.
         */
        Dataset2D(std::vector<double> x, std::vector<double> y, std::vector<double> xerr, std::vector<double> yerr) noexcept;

        Dataset2D(const SimpleDataset& data);

        /**
         * @brief Construct a new dataset from an input file.
         */
        Dataset2D(std::string path);

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

        // Get the fourth column.
        [[nodiscard]] const ConstColumn<double> xerr() const {return ConstColumn<double>(data, N, M, 3);}

        // Get the fourth column.
        [[nodiscard]] Column<double> xerr() {return Column<double>(data, N, M, 3);}

        // Get the ith value in the fourth column.
        [[nodiscard]] const double& xerr(unsigned int i) const {return index(i, 3);}

        // Get the ith value in the fourth column.
        [[nodiscard]] double& xerr(unsigned int i) {return index(i, 3);}
};

// Object conversion between Dataset2D and SimpleDataset is often used. This check ensures that the conversion is safe.
static_assert(sizeof(Dataset2D) == sizeof(SimpleDataset), "Object conversion is broken.");