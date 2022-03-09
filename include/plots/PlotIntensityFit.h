#pragma once

#include "plots/Plot.h"
#include "ScatteringHistogram.h"
#include "fitter/SimpleIntensityFitter.h"
#include "settings.h"

#include <memory.h>
#include <string.h>
#include <vector>

#include <TLegend.h>
#include <TH1D.h>
#include <TLine.h>

using std::unique_ptr, std::shared_ptr, std::string, std::vector;

namespace plots {
  class PlotIntensityFit : public Plot {
    public:
      PlotIntensityFit(SimpleIntensityFitter& fitter) : Plot(), fitter(fitter) {}

      ~PlotIntensityFit() override = default;

      void save(std::string path) const override;

    private:
      SimpleIntensityFitter& fitter;
  };
}