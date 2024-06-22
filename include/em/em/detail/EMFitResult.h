#pragma once

#include <fitter/FitResult.h>

namespace fitter {
    class EMFitResult : public FitResult {
        public:
            using FitResult::FitResult;
            ~EMFitResult() override = default;

            [[nodiscard]] std::string to_string() const noexcept override;

            double level = 0;
            double mass = 0;

            struct EMFitPlots {SimpleDataset chi2_full, chi2_limited, chi2_minimum, mass_full, mass_limited, mass_minimum;} em_plots;
    };
}