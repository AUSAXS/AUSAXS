#pragma once

#include <fitter/FitterFwd.h>

#include <mini/detail/Landscape.h>
#include <mini/detail/Result.h>
#include <dataset/SimpleDataset.h>
#include <utility/observer_ptr.h>

#include <string>

namespace fitter {
    class FitResult : public mini::Result {
        public:
            FitResult() noexcept = default;
            FitResult(const mini::Result& res, double chi2, unsigned int dof) noexcept;
            ~FitResult() override = default;
            
            /**
             * @brief Add the parameters from another fit to this one. Each parameter will count as an additional degree of freedom. 
             *
             * @param fit The fit to add.
             * @param front If true, the parameters will be added to the front of the parameter list. Otherwise, they will be added to the back.
             */
            void add_fit(observer_ptr<FitResult> fit, bool front = false) noexcept;

            /**
             * @brief Get a string representation of this object. 
             */
            [[nodiscard]] virtual std::string to_string() const noexcept;

            Dataset curves; // | q | data | interpolated model | residuals |
            mini::Landscape evaluated_points;
            unsigned int dof;
    };
}