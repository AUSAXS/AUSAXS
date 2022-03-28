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
         * Create a new dataset based on a list of x and y coordinates. 
         */
        Dataset(const std::vector<double>& x, const std::vector<double>& y);

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
        Dataset& reduce(unsigned int target);

        /**
         * @brief Impose limits on the data. All points with an x-value outside this range will be removed. 
         * 
         * @param limits The new limits. 
         * @return The modified dataset. 
         */
        Dataset& limit(const Limit& limits);

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

        std::unique_ptr<TGraph> plot() const;

        std::string xlabel = "x";
        std::string ylabel = "y";
        std::vector<double> x; // The x coordinates.
        std::vector<double> y; // The y coordinates.

    private:
        void check_sizes() const;
};