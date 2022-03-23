#pragma once

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
      PlotIntensity(const ScatteringHistogram& d);

      /**
       * @brief Move constructor.
       * 
       * @param d The ScatteringHistogram to be plotted. 
       */
      PlotIntensity(ScatteringHistogram&& d) : Plot(), d(std::move(d)) {}

      /**
       * @brief Plot an additional data set as points. 
       */
      void plot_intensity(const Dataset& data, EColor color = kBlack);

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

      void plot_intensity();
  };
}