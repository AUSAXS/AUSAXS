#pragma once

#include "plots/Plot.h"
#include "ScatteringHistogram.h"
#include "settings.h"

namespace plots {
  class PlotIntensity : public Plot {
    public:
      /**
       * @brief Copy constructor.
       * 
       * @param d The ScatteringHistogram to be plotted. 
       */
      PlotIntensity(const ScatteringHistogram& d) : Plot(), d(d) {}

      /**
       * @brief Move constructor.
       * 
       * @param d The ScatteringHistogram to be plotted. 
       */
      PlotIntensity(ScatteringHistogram&& d) : Plot(), d(std::move(d)) {}

      /**
       * @brief Destructor.
       */
      ~PlotIntensity() override = default;

      void save(std::string path) const override;

    private:
      const ScatteringHistogram d;
  };
}