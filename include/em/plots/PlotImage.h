#pragma once

#include <plots/Plot.h>

namespace ausaxs::em {class Image;}
namespace ausaxs::plots {
    /**
     * @brief Plot a specific \class Image object. 
     */
    class PlotImage : public Plot {
        public:
            /**
             * @brief Constructor.
             * 
             * @param image The image to be plotted. 
             */
            PlotImage(const em::Image& image, const PlotOptions& options);

            /**
             * @brief Destructor. 
             */
            ~PlotImage() override;

            /**
             * @brief Overlay the atomic locations on top of the plot. 
             * 
             * @param cutoff The cutoff value of the electron density. 
             */
            PlotImage& plot_atoms(const em::Image& image, double cutoff);

			/**
			 * @brief Plot and save the input Image at the specified location. 
			 * 	      This is a convenient shortcut for quickly creating a plot of a single Image. 
			 */
            static void quick_plot(const em::Image& image, const PlotOptions& options, const io::File& path);

        private:
            void plot(const em::Image& image, const PlotOptions& options);
    };
}