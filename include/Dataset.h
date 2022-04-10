#pragma once

#include <vector>
#include <string>
#include <memory>

#include <TGraph.h>

#include <data/Axis.h>

/**
 * @brief A representation of a set of 2D data. 
 */
class Dataset {
    public: 
        /**
         * @brief Default constructor.
         */
        Dataset() {}

        /**
         * @brief Constructor.
         * 
         * Create a new dataset based on an input file. Format is assumed to be x | y | yerr
         * 
         * @param file 
         */
        Dataset(const std::string file);

        /**
         * @brief Constructor. 
         * 
         * Create a new dataset based on a list of x and y coordinates. 
         */
        Dataset(const std::vector<double>& x, const std::vector<double>& y);

        /**
         * @brief Constructor. 
         * 
         * Create a new dataset based on a list of x and y coordinates, along with an error on the latter. 
         */
        Dataset(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& yerr);

        /**
         * @brief Constructor. 
         * 
         * Create a new dataset based on a list of x and y coordinates, along with errors for both.
         */
        Dataset(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& yerr, const std::vector<double>& xerr);

        /**
         * @brief Constructor. 
         * 
         * Create a new labelled dataset based on a list of x and y coordinates. 
         */
        Dataset(const std::vector<double>& x, const std::vector<double>& y, const std::string xlabel, const std::string ylabel);

        /**
         * @brief Destructor. 
         */
        ~Dataset() = default;

        /**
         * @brief Reduce the number of data points to the specified amount. 
         * 
         * @return The modified dataset. 
         */
        void reduce(unsigned int target, bool log = false);

        /**
         * @brief Impose limits on the data. All points with an x-value outside this range will be removed. 
         * 
         * @param limits The new limits. 
         * @return The modified dataset. 
         */
        void limit(const Limit& limits);

        /**
         * @brief Get the number of data points. 
         */
        std::size_t size() const;

        /**
         * @brief Get a set of data based on its label. 
         * This also serves as a consistency check. 
         * 
         * @param label The name of the data. 
         */
        std::vector<double>& get(const std::string label);

        /**
         * @brief Scale all errors by some common factor. 
         */
        void scale_errors(double factor);

        /**
         * @brief Scale the y-values (and their associated errors) by some common factor.
         */
        void scale_y(double factor);

        std::unique_ptr<TGraph> plot() const;

        std::string xlabel = "x";
        std::string ylabel = "y";
        std::string xerrlabel = "xerr";
        std::string yerrlabel = "yerr";
        std::vector<double> x;    // The x coordinates.
        std::vector<double> y;    // The y coordinates.
        std::vector<double> xerr; // The error in the x coordinates
        std::vector<double> yerr; // The error in the y coordinates
        bool draw_as_line = true;

    private:
        void validate_sizes() const;

        void read(const std::string file);
};

class SAXSDataset : public Dataset {
    using Dataset::Dataset; // inherit constructors

    public:
        /**
         * @brief Generate errors for the y-values mimicking what one would find experimentally. 
         */
        void simulate_errors();

        /**
         * @brief Set the resolution of this dataset. 
         */
        void set_resolution(unsigned int resolution);

    private:
        unsigned int resolution = 0;
};