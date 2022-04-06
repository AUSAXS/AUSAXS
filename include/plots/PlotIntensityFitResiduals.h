#pragma once

#include <plots/Plot.h>
#include <ScatteringHistogram.h>
#include <fitter/IntensityFitter.h>
#include <settings.h>

#include <memory.h>
#include <string.h>
#include <vector>

namespace plots {
  /**
   * @brief \class PlotIntensityFitResiduals
   * 
   * Plot the residuals of the fitted scattering curve. 
   * Remember to set the correct ScatteringPlot with the optimized values in the fitter before using this class. 
   */
  class PlotIntensityFitResiduals : public Plot {
    public:

      /**
       * @brief Constructor.
       * 
       * @param fitter The fit to plot. Remember to update it with the optimized values before creating an instance of this class. 
       */
      PlotIntensityFitResiduals(SimpleIntensityFitter& fitter);

      /**
       * @brief Constructor.
       * 
       * @param fitter The fit to plot. Remember to update it with the optimized values before creating an instance of this class. 
       */
      PlotIntensityFitResiduals(const Fitter::Fit& fit);

      /**
       * @brief Constructor.
       * 
       * @param fitter The fit to plot. Remember to update it with the optimized values before creating an instance of this class. 
       */
      PlotIntensityFitResiduals(const std::shared_ptr<Fitter::Fit> fit);

      /**
       * @brief Destructor.
       */
      ~PlotIntensityFitResiduals() override = default;

      /**
       * @brief Create and save the plot at the given path. 
       * 
       * @param path Save location and format. 
       */
      void save(std::string path) const override;

    private:
      std::unique_ptr<TCanvas> canvas;

      void plot(const std::shared_ptr<TGraph> graph) const;

      void prepare_canvas();
  };
}