#pragma once

#include <TCanvas.h>
#include <TPad.h>
#include <TH2D.h>

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
            ~PlotImage() = default;

            /**
             * @brief Save this image at the given location in the specified format. 
             * 
             * @param path The path & format of the image. 
             */
            void save(std::string path) const override;

            /**
             * @brief Overlay the atomic locations on top of the plot. 
             * 
             * @param cutoff The cutoff value of the electron density. 
             */
            void plot_atoms(double cutoff = 0.1) const;

        private:
            std::unique_ptr<TCanvas> canvas;
            std::unique_ptr<TPad> pad1;
            std::unique_ptr<TPad> pad2;
            const em::Image& image;

            std::unique_ptr<TH2D> plot_hist() const;
    };
}