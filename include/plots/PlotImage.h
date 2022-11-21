#pragma once

#include <plots/Plot.h>
#include <em/Image.h>

namespace plots {
    /**
     * @brief \class PlotImage
     * 
     * Plot a specific \class Image object. 
     */
    class PlotImage : public Plot {
        public:
            /**
             * @brief Constructor.
             * 
             * @param image The image to be plotted. 
             */
            PlotImage(const em::Image& image);

            /**
             * @brief Destructor. 
             */
            ~PlotImage() override;

            /**
             * @brief Overlay the atomic locations on top of the plot. 
             * 
             * @param cutoff The cutoff value of the electron density. 
             */
            void plot_atoms(const em::Image& image, double cutoff);

			/**
			 * @brief Plot and save the input Image at the specified location. 
			 * 	      This is a convenient shortcut for quickly creating a plot of a single Image. 
			 */
            static void quick_plot(const em::Image& image, std::string path);

        private:
            void plot(const em::Image& image);
    };
}