#pragma once

#include "plots/Plot.h"
#include "ScatteringHistogram.h"
#include "fitter/IntensityFitter.h"
#include "settings.h"

#include <memory.h>
#include <string.h>
#include <vector>

namespace plots {
  class PlotIntensityFitResiduals : public Plot {
    public:
      PlotIntensityFitResiduals(SimpleIntensityFitter& fitter) : Plot(), fitter(fitter) {}
      ~PlotIntensityFitResiduals() override = default;

      void save(std::string path) const override;

    private:
      SimpleIntensityFitter& fitter;
  };
}