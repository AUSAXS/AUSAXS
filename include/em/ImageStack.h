#pragma once

#include <vector>
#include <string>

#include <em/CullingStrategy.h>
#include <em/Datatypes.h>
#include <em/Image.h>
#include <data/Protein.h>
#include <hydrate/Grid.h>
#include <histogram/ScatteringHistogram.h>
#include <fitter/SimpleIntensityFitter.h>
#include <fitter/IntensityFitter.h>
#include <utility/Multiset.h>
#include <minimizer/Minimizer.h>

namespace em {
    class PartialHistogramManager;

    /**
     * @brief \class ImageStack
     * 
     * A representation of a stack of images. 
     */
    class ImageStack {
        public:
            class EMFit : public Fit {using Fit::Fit;};

            struct Landscape {
                Landscape() {}
                Landscape(const EMFit& fit, const Dataset2D& contour) : fit(fit), contour(contour) {}
                Landscape(Dataset2D&& contour) : contour(std::move(contour)) {}

                EMFit fit;
                Dataset2D contour;
            };

            /**
             * @brief Constructor.
             * 
             * @param file Path to the input EM data file. 
             */
            ImageStack(std::string file, unsigned int resolution = 0);

            /**
             * @brief Constructor.
             * 
             * @param resolution 
             */
            ImageStack(const std::vector<Image>& images, unsigned int resolution = 0);

            /**
             * @brief Destructor.
             */
            ~ImageStack();;

            /**
             * @brief Fit the cutoff value with the input experimental data file. 
             * 
             * @param file Path to the measurement file. 
             * @param param The cutoff parameter.
             */
            std::shared_ptr<EMFit> fit(std::string file, mini::Parameter param);

            /**
             * @brief Fit the cutoff value with the input experimental data file. 
             * 
             * @param file Path to the measurement file. 
             */
            std::shared_ptr<EMFit> fit(std::string file);

            /**
             * @brief Fit the cutoff value with the input histogram. 
             * 
             * @param h The histogram to fit to.  
             * @param param The cutoff parameter.
             */
            std::shared_ptr<EMFit> fit(const hist::ScatteringHistogram& h, mini::Parameter param);

            /**
             * @brief Fit the cutoff value with the input histogram. 
             * 
             * @param h The histogram to fit to.  
             */
            std::shared_ptr<EMFit> fit(const hist::ScatteringHistogram& h);

            /**
             * @brief Perform a scan of the cutoff values. 
             * 
             * @param points The cutoff values to be evaluated. 
             * @param h The histogram to be fitted. 
             * 
             * @return A Dataset containing the scanned cutoff values and their corresponding chi2 values. 
             */
            Landscape cutoff_scan(const Axis& points, const hist::ScatteringHistogram& h);

            /**
             * @brief Perform a scan of the cutoff values. 
             * 
             * @param points The number of points.
             * @param h The histogram to be fitted. 
             * 
             * @return A Dataset containing the scanned cutoff values and their corresponding chi2 values. 
             */
            Landscape cutoff_scan(unsigned int points, const hist::ScatteringHistogram& h);

            /**
             * @brief Perform a scan of the cutoff values. 
             * 
             * @param points The cutoff values to be evaluated. 
             * @param h The measurement file to compare with. 
             * 
             * @return A Dataset containing the scanned cutoff values and their corresponding chi2 values. 
             */
            Landscape cutoff_scan(const Axis& points, std::string file);

            /**
             * @brief Perform a scan of the cutoff values. 
             * 
             * @param points The number of points.
             * @param h The measurement file to compare with. 
             * 
             * @return A Dataset containing the scanned cutoff values and their corresponding chi2 values. 
             */
            Landscape cutoff_scan(unsigned int points, std::string file);

            /**
             * @brief Perform a scan & fit of the cutoff values. 
             * 
             * @param points The cutoff values to be evaluated. 
             * @param h The histogram to be fitted. 
             * 
             * @return A Landscape containing both the fit and scan.
             */
            Landscape cutoff_scan_fit(const Axis& points, const hist::ScatteringHistogram& h);

            /**
             * @brief Perform a scan & fit of the cutoff values. 
             * 
             * @param points The number of points.
             * @param h The histogram to be fitted. 
             * 
             * @return A Landscape containing both the fit and scan.
             */
            Landscape cutoff_scan_fit(unsigned int points, const hist::ScatteringHistogram& h);

            /**
             * @brief Perform a scan & fit of the cutoff values. 
             * 
             * @param points The number of points.
             * @param file The file to be fitted.
             * 
             * @return A Landscape containing both the fit and scan.
             */
            Landscape cutoff_scan_fit(unsigned int points, std::string file);

            /**
             * @brief Perform a scan & fit of the cutoff values. 
             * 
             * @param points The cutoff values to be evaluated. 
             * @param file The file to be fitted.
             * 
             * @return A Landscape containing both the fit and scan.
             */
            Landscape cutoff_scan_fit(const Axis& points, std::string file);

            /**
             * @brief Get a specific Image stored in this object. 
             * 
             * @param layer The vertical location of the Image. 
             */
            Image& image(unsigned int layer);

            /**
             * @brief Get a specific Image stored in this object. 
             * 
             * @param layer The vertical location of the Image. 
             */
            const Image& image(unsigned int layer) const;

            /**
             * @brief Prepare a ScatteringHistogram based on this object. 
             */
            hist::ScatteringHistogram get_histogram(double cutoff) const;

            /**
             * @brief Count the number of voxels for a given cutoff.
             */
            unsigned int count_voxels(double cutoff) const;

            /**
             * @brief Get the fitted ScatteringHistogram.
             */
            hist::ScatteringHistogram get_histogram(const std::shared_ptr<EMFit> res) const;

            /**
             * @brief Get the protein generated with the chosen cutoff value.
             */
            std::shared_ptr<Protein> get_protein(double cutoff) const;

            /**
             * @brief Create a new Grid based on this object. 
             * 
             * @param cutoff The cutoff value. If positive, atoms will be generated at all pixel values higher than this. If negative, they will be generated at pixels lower than this. 
             */
            std::unique_ptr<Grid> create_grid(double cutoff) const;

            /**
             * @brief Get the header of the input file. 
             */
            std::shared_ptr<ccp4::Header> get_header() const;

            /**
             * @brief Get the number of images stored in this object.
             */
            size_t size() const;

            /**
             * @brief Get a reference to all images stored in this object. 
             */
            const std::vector<Image>& images() const;

            /**
             * @brief Save this structure as a .pdb file. 
             * 
             * @param path Path to save location.
             * @param cutoff The cutoff value. If positive, atoms will be generated at all pixel values higher than this. If negative, they will be generated at pixels lower than this. 
             */
            void save(std::string path, double cutoff) const;

            /**
             * @brief Get the limits on the q-values generated by discretizing the model scattering curve.
             * 
             * More specifically the upper limit will be $q_{max} = \frac{2\pi}{d}$ where $d$ is the resolution assumed to be the only number present in the filename. 
             */
            Limit get_limits() const;

            /**
             * @brief Get the mean density.
             */
            double mean() const;

            ObjectBounds3D minimum_volume(double cutoff);

            /**
             * @brief Get the cutoff corresponding to a PyMOL level. This is just the number of sigmas of the root-mean-square deviation.
             */
            double level(double sigma) const;

            /**
             * @brief Calculate the root-mean-square of this map. 
             */
            double rms() const;

            /**
             * @brief Get the histogram manager.
             */
            std::shared_ptr<em::PartialHistogramManager> get_histogram_manager() const;

            /**
             * @brief Get the fitted water scaling factors.
             */
            const std::vector<mini::FittedParameter>& get_fitted_water_factors() const;

            /**
             * @brief Get the fitted water scaling factors as a dataset.
             */
            SimpleDataset get_fitted_water_factors_dataset() const;

        private:
            std::string filename;
            std::shared_ptr<ccp4::Header> header;               // The header of the map file. 
            std::vector<Image> data;                            // The actual image data. 
            unsigned int resolution;                            // The resolution of the map //? Deprecated?
            unsigned int size_x, size_y, size_z;                // The number of pixels in each dimension.
            std::shared_ptr<em::PartialHistogramManager> phm;   // The histogram manager. Manages both the backing protein & its scattering curve. 
            std::vector<mini::FittedParameter> water_factors;   // If hydration is enabled, the fitted water scaling factors will be recorded here.

            /**
             * @brief Update the cutoff sections that will be used.
             */
            void update_charge_levels(Limit limit) const noexcept;

            /**
             * @brief Determines the minimum bounds necessariy to describe the map for the given cutoff.
             * 
             * @param min_val The smallest possible value. Must be positive - sign is determined based on the determine_staining() method.
             * TODO Determine a better, more dynamic approach to determining this minimum. 
             */
            void determine_minimum_bounds(double min_val);

            /**
             * @brief Check if the extension is valid.
             *        Throws an exception if not.
             */
            void validate_extension(string file) const;
            
            void read(std::ifstream& istream, size_t byte_size);

            /**
             * @brief Get the data byte size of the CCP file. 
             */
            size_t get_byte_size() const;

            /**
             * @brief A helper function for the fitting methods. This performs the actual fit. 
             * 
             * @param fitter The fitter object to fit. 
             */
            std::shared_ptr<EMFit> fit_helper(std::shared_ptr<SimpleIntensityFitter> fitter, mini::Parameter param = {});

            /**
             * @brief A helper function for the cutoff scanning method.
             * 
             * @param points The range to scan.
             * @param fitter The fitting object.
             */
            Landscape cutoff_scan_helper(const Axis& points, std::shared_ptr<SimpleIntensityFitter> fitter);

            /**
             * @brief A helper function for the cutoff scan & fit method.
             * 
             * @param points The range to scan.
             * @param fitter The fitting object.
             */
            Landscape cutoff_scan_fit_helper(const Axis& points, std::shared_ptr<SimpleIntensityFitter> fitter);

            /**
             * @brief Prepare the fitting function. 
             *        Note that the lifetime of the returned function is the same as that of the fitter.
             */
            std::function<double(const double*)> prepare_function(std::shared_ptr<SimpleIntensityFitter> fitter);

            float& index(unsigned int x, unsigned int y, unsigned int z);

            float index(unsigned int x, unsigned int y, unsigned int z) const;
    };
}