#pragma once

#include <plots/Plot.h>
#include <ScatteringHistogram.h>
#include <fitter/SimpleIntensityFitter.h>
#include <settings.h>

#include <memory.h>
#include <string.h>
#include <vector>

namespace plots {

  /**
   * @brief \class PlotIntensityFit
   * 
   * Plot both the measured and fitted scattering curve. 
   * Remember to set the correct ScatteringPlot with the optimized values in the fitter before using this class. 
   */
  class PlotIntensityFit : public Plot {
    public:
      /**
       * @brief Constructor.
       * 
       * @param fitter The fit to plot. Remember to update it with the optimized values before creating an instance of this class. 
       */
      PlotIntensityFit(SimpleIntensityFitter& fitter);

      /**
       * @brief Constructor.
       * 
       * @param fitter The fit to plot. Remember to update it with the optimized values before creating an instance of this class. 
       */
      PlotIntensityFit(const Fitter::Fit& fit);

      /**
       * @brief Constructor.
       * 
       * @param fitter The fit to plot. Remember to update it with the optimized values before creating an instance of this class. 
       */
      PlotIntensityFit(const std::shared_ptr<Fitter::Fit> fit);

      /**
       * @brief Destructor.
       */
      ~PlotIntensityFit() override = default;

      /**
       * @brief Create and save the plot at the given path. 
       * 
       * @param path Save location and format. 
       */
      void save(std::string path) const override;

    private:
      std::unique_ptr<TCanvas> canvas;

      void plot(const std::vector<std::shared_ptr<TGraph>>& graphs) const;

      void prepare_canvas();
  };
}