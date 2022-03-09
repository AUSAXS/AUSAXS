#pragma once

#include <TCanvas.h>
#include <TPad.h>
#include <TH2D.h>

#include <plots/Plot.h>
#include <em/Image.h>

namespace plots {
    class PlotImage : public Plot {
        public:
            PlotImage(const em::Image& image);

            ~PlotImage() = default;

            void save(std::string path) const override;

            void plot_atoms(double cutoff = 0.1) const;

        private:
            std::unique_ptr<TCanvas> canvas;
            std::unique_ptr<TPad> pad1;
            std::unique_ptr<TPad> pad2;
            const em::Image& image;

            std::unique_ptr<TH2D> plot_hist() const;
    };
}