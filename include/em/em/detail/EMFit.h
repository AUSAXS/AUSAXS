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
    };
}