#pragma once

#include <em/detail/ImageStackBase.h>

#include <io/IOFwd.h>
#include <hist/HistFwd.h>
#include <mini/MiniFwd.h>
#include <fitter/FitterFwd.h>
#include <dataset/DatasetFwd.h>
#include <em/detail/EMInternalFwd.h>

#include <utility/Observable.h>

#include <functional>

namespace em {
    /**
     * @brief Extends the ImageStackBase class with fitting functionalities. 
     */
    class ImageStack : public ImageStackBase {
        public:            
            ImageStack() = default;

            /**
             * @brief Create a new ImageStack from an input EM data file.
             * 
             * @param file Path to the input EM data file. 
             */
            ImageStack(const io::ExistingFile& file);

            /**
             * @brief Create a new ImageStack from a list of images.
             * 
             * @param images The images for this stack.
             */
            ImageStack(const std::vector<Image>& images);

            ~ImageStack() override;

            /**
             * @brief Get the mass corresponding to the given cutoff level. 
             * 
             * Complexity: O(n) where n is the number of atoms in the stack, with a significant prefactor. 
             * 
             * @return The mass in kDa.
             */
            double get_mass(double cutoff) const; 

            /**
             * @brief Fit the cutoff value with the input experimental data file. 
             * 
             * @param file Path to the measurement file. 
             * @param param The cutoff parameter.
             */
            std::unique_ptr<fitter::EMFitResult> fit(const io::ExistingFile& file, mini::Parameter& param);

            /**
             * @brief Fit the cutoff value with the input experimental data file. 
             * 
             * @param file Path to the measurement file. 
             */
            std::unique_ptr<fitter::EMFitResult> fit(const io::ExistingFile& file);

            // /**
            //  * @brief Fit the cutoff value with the input histogram. 
            //  * 
            //  * @param h The histogram to fit to.  
            //  * @param param The cutoff parameter.
            //  */
            // std::unique_ptr<fitter::EMFitResult> fit(std::unique_ptr<hist::ICompositeDistanceHistogram> h, mini::Parameter& param);

            // /**
            //  * @brief Fit the cutoff value with the input histogram. 
            //  * 
            //  * @param h The histogram to fit to.  
            //  */
            // std::unique_ptr<fitter::EMFitResult> fit(std::unique_ptr<hist::ICompositeDistanceHistogram> h);

            // /**
            //  * @brief Perform a scan of the cutoff values. 
            //  * 
            //  * @param points The cutoff values to be evaluated. 
            //  * @param h The histogram to be fitted. 
            //  * 
            //  * @return A Dataset containing the scanned cutoff values and their corresponding chi2 values. 
            //  */
            // mini::Landscape cutoff_scan(const Axis& points, std::unique_ptr<hist::ICompositeDistanceHistogram> h);

            // /**
            //  * @brief Perform a scan of the cutoff values. 
            //  * 
            //  * @param points The number of points.
            //  * @param h The histogram to be fitted. 
            //  * 
            //  * @return A Dataset containing the scanned cutoff values and their corresponding chi2 values. 
            //  */
            // mini::Landscape cutoff_scan(unsigned int points, std::unique_ptr<hist::ICompositeDistanceHistogram> h);

            // /**
            //  * @brief Perform a scan of the cutoff values. 
            //  * 
            //  * @param points The cutoff values to be evaluated. 
            //  * @param h The measurement file to compare with. 
            //  * 
            //  * @return A Dataset containing the scanned cutoff values and their corresponding chi2 values. 
            //  */
            // mini::Landscape cutoff_scan(const Axis& points, const io::ExistingFile& file);

            // /**
            //  * @brief Perform a scan of the cutoff values. 
            //  * 
            //  * @param points The number of points.
            //  * @param h The measurement file to compare with. 
            //  * 
            //  * @return A Dataset containing the scanned cutoff values and their corresponding chi2 values. 
            //  */
            // mini::Landscape cutoff_scan(unsigned int points, const io::ExistingFile& file);

            // /**
            //  * @brief Perform a scan & fit of the cutoff values. 
            //  * 
            //  * @param points The cutoff values to be evaluated. 
            //  * @param h The histogram to be fitted. 
            //  * 
            //  * @return A Landscape containing both the fit and scan.
            //  */
            // std::pair<fitter::EMFitResult, mini::Landscape> cutoff_scan_fit(const Axis& points, std::unique_ptr<hist::ICompositeDistanceHistogram> h);

            // /**
            //  * @brief Perform a scan & fit of the cutoff values. 
            //  * 
            //  * @param points The number of points.
            //  * @param h The histogram to be fitted. 
            //  * 
            //  * @return A Landscape containing both the fit and scan.
            //  */
            // std::pair<fitter::EMFitResult, mini::Landscape> cutoff_scan_fit(unsigned int points, std::unique_ptr<hist::ICompositeDistanceHistogram> h);

            // /**
            //  * @brief Perform a scan & fit of the cutoff values. 
            //  * 
            //  * @param points The number of points.
            //  * @param file The file to be fitted.
            //  * 
            //  * @return A Landscape containing both the fit and scan.
            //  */
            // std::pair<fitter::EMFitResult, mini::Landscape> cutoff_scan_fit(unsigned int points, const io::ExistingFile& file);

            // /**
            //  * @brief Perform a scan & fit of the cutoff values. 
            //  * 
            //  * @param points The cutoff values to be evaluated. 
            //  * @param file The file to be fitted.
            //  * 
            //  * @return A Landscape containing both the fit and scan.
            //  */
            // std::pair<fitter::EMFitResult, mini::Landscape> cutoff_scan_fit(const Axis& points, const io::ExistingFile& file);

            /**
             * @brief Get the fitted water scaling factors.
             */
            const std::vector<mini::FittedParameter>& get_fitted_water_factors() const;

            /**
             * @brief Get the fitted water scaling factors as a dataset.
             */
            SimpleDataset get_fitted_water_factors_dataset() const;

            /**
             * @brief Get a new progress observer.
             * 
             * This will be notified of the progress of the fitting process in the form of a iteration counter.
             * Typically this will be in the range [0, 2*max_iterations], though it can be higher for some datasets.
             */
            auto get_progress_observer() {return progress.make_observer();}

        private: 
            std::vector<mini::FittedParameter> water_factors;   // If hydration is enabled, the fitted water scaling factors will be recorded here.
            std::vector<detail::ExtendedLandscape> evals;       // The evaluated points.
            utility::Observable<int> progress;                  // The progress of the fitting process.

            /**
             * @brief Update the cutoff sections that will be used.
             */
            void update_charge_levels(const Limit& limit) const noexcept;

            /**
             * @brief Prepare the fitting function. 
             *        Note that the lifetime of the returned function is the same as that of the fitter.
             */
            std::function<double(std::vector<double>)> prepare_function(std::shared_ptr<fitter::SmartFitter> fitter);

            // /**
            //  * @brief A helper function for the fitting methods. This performs the actual fit. 
            //  * 
            //  * @param fitter The fitter object to fit. 
            //  */
            // std::unique_ptr<fitter::EMFitResult> fit_helper(std::shared_ptr<fitter::SmartFitter> fitter);

            /**
             * @brief A helper function for the fitting methods. This performs the actual fit. 
             * 
             * @param fitter The fitter object to fit. 
             */
            std::unique_ptr<fitter::EMFitResult> fit_helper(std::shared_ptr<fitter::SmartFitter> fitter, mini::Parameter& param);

            // /**
            //  * @brief A helper function for the cutoff scanning method.
            //  * 
            //  * @param points The range to scan.
            //  * @param fitter The fitting object.
            //  */
            // mini::Landscape cutoff_scan_helper(const Axis& points, std::shared_ptr<fitter::SmartFitter> fitter);

            // /**
            //  * @brief A helper function for the cutoff scan & fit method.
            //  * 
            //  * @param points The range to scan.
            //  * @param fitter The fitting object.
            //  */
            // std::pair<fitter::EMFitResult, mini::Landscape> cutoff_scan_fit_helper(const Axis& points, std::shared_ptr<fitter::SmartFitter> fitter);
    };
}