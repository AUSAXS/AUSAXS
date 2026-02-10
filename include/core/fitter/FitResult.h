// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <fitter/FitterFwd.h>

#include <mini/detail/Landscape.h>
#include <mini/detail/Result.h>
#include <dataset/NamedDataset.h>
#include <utility/observer_ptr.h>

#include <string>

namespace ausaxs::fitter {
    class FitResult : public mini::Result {
        public:
            FitResult() noexcept = default;
            FitResult(const mini::Result& res, unsigned int dof) noexcept;
            ~FitResult() override = default;
            
            /**
             * @brief Add the parameters from another fit to this one. Each parameter will count as an additional degree of freedom. 
             *
             * @param fit The fit to add.
             * @param front If true, the parameters will be added to the front of the parameter list. Otherwise, they will be added to the back.
             */
            void add_fit(observer_ptr<FitResult> fit, bool front = false) noexcept;

            /**
             * @brief Set the data curves for this fit. 
             *
             * @throw std::invalid_argument If the number of columns is not 5 or the column names are not as expected. 
             */
            void set_data_curves(NamedDataset&& curves);

            /**
             * @brief Set the data curves for this fit. 
             */
            void set_data_curves(std::vector<double>&& q, std::vector<double>&& data, std::vector<double>&& data_err, std::vector<double>&& model, std::vector<double>&& residuals);

            /**
             * @brief Get a string representation of this object. 
             */
            [[nodiscard]] virtual std::string to_string() const noexcept;

            NamedDataset curves; // | q | data | data_err | interpolated model | residuals |
            mini::Landscape evaluated_points;
            unsigned int dof;
    };
}