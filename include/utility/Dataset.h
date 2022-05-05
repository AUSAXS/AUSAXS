#pragma once

#include <vector>
#include <string>
#include <memory>
#include <any>

#include <TGraph.h>

#include <data/Axis.h>
#include <plots/PlotOptions.h>

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
         * Create a new empty dataset with the given labels.
         */
        Dataset(std::string xlabel, std::string ylabel);

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
        Dataset(const std::vector<double>& x, const std::vector<double>& y, std::string xlabel, std::string ylabel);

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

        /**
         * @brief Set the normalization of the y-values. The first y-value will be fixed to this. 
         */
        void normalize(double y0);

        /**
         * @brief Plot this dataset.
         */
        std::unique_ptr<TGraph> plot() const;

        /**
         * @brief Set the plot options of this dataset.
         */
        void set_plot_options(const plots::PlotOptions& options);

        /**
         * @brief Set the plot options of this dataset.
         */
        void set_plot_options(const std::map<std::string, std::any>& options);

        /**
         * @brief Set the plot options of this dataset.
         * 
         * @param style The plotting style. Should be one of the accepted variations of "markers", "errors", or "line". 
         * @param options The other plot options.
         */
        void set_plot_options(std::string style, std::map<std::string, std::any> options = {});

        /**
         * @brief Draw a dataset with its currently set plot options.
         */
        static void draw(const Dataset& data);

        /**
         * @brief Write this dataset to the specified file. 
         */
        void save(std::string path) const;

        std::string xlabel = "x";
        std::string ylabel = "y";
        std::string xerrlabel = "xerr";
        std::string yerrlabel = "yerr";
        std::vector<double> x;    // The x coordinates.
        std::vector<double> y;    // The y coordinates.
        std::vector<double> xerr; // The error in the x coordinates
        std::vector<double> yerr; // The error in the y coordinates
        plots::PlotOptions plot_options; 

    private:
        /**
         * @brief Validate the sizes of the stored vectors. 
         */
        void validate_sizes() const;

        /**
         * @brief Read a dataset from a file.
         * 
         * @param file Path to the input file. 
         */
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
         * @brief Simulate Gaussian noise on the y-values based on the errors. 
         */
        void simulate_noise();

        /**
         * @brief Set the resolution of this dataset. 
         */
        void set_resolution(unsigned int resolution);

    private:
        unsigned int resolution = 0;
};