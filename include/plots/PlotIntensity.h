#pragma once

#include <fitter/Fitter.h>
#include <plots/Plot.h>
#include <ScatteringHistogram.h>
#include <settings.h>

#include <TH2D.h>

namespace plots {
  class PlotIntensity : public Plot {
    public:
      /**
       * @brief Copy constructor.
       * 
       * @param d The ScatteringHistogram to be plotted. 
       */
      PlotIntensity(const ScatteringHistogram& d, int color = kBlack);

      /**
       * @brief Move constructor.
       * 
       * @param d The ScatteringHistogram to be plotted. 
       */
      PlotIntensity(ScatteringHistogram&& d, int color = kBlack);

      /**
       * @brief Constructor.
       * 
       * @param d The dataset to be plotted.
       */
      PlotIntensity(const SAXSDataset& d, int color = kBlack);

      /**
       * @brief Plot an additional data set as points. 
       */
      void plot_intensity(const Dataset& data, int color = kBlack, double alpha = 1);

      /**
       * @brief Plot the result of a fit. 
       */
      void plot_intensity(const std::shared_ptr<Fitter::Fit> fit, int color = kBlack, double alpha = 1);

      /**
       * @brief Destructor.
       */
      ~PlotIntensity() override = default;

      void plot_guinier_approx();

      void save(std::string path) const override;

    private:
      const ScatteringHistogram d;
      unique_ptr<TCanvas> canvas;
      unique_ptr<TPad> linpad;
      unique_ptr<TPad> logpad;
      double ymin, ymax;

      void initial_intensity_plot(int color);

      void initial_intensity_plot(const Dataset& data, int color);

      void prepare_canvas();
  };
}