// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <fitter/detail/LinearLeastSquares.h>
#include <dataset/DatasetFwd.h>
#include <hist/HistFwd.h>

#include <vector>

namespace ausaxs::fitter {
    /**
     * @brief A simple linear least-squares fitter for fitting the linear relationship y = ax+b.
     */
    class LinearFitter : public detail::LinearLeastSquares {
        public:
            virtual ~LinearFitter() override = default;
            LinearFitter(LinearFitter&&) noexcept;
            LinearFitter& operator=(LinearFitter&&) noexcept;

            /**
             * @brief Prepare a fit of the measured values in @a input to a model to be defined later. 
             */
            LinearFitter(const SimpleDataset& data);

            /**
             * @brief Prepare a fit of the measured values in @a input to the model described by @a h.
             */
            LinearFitter(const SimpleDataset& data, std::unique_ptr<hist::DistanceHistogram> model);

            /**
             * @brief Get the dataset being fitted. 
             */
            SimpleDataset get_data() const;

			/**
			 * @brief Set the scattering histogram used for the fit. 
			 */
			virtual void set_model(std::unique_ptr<hist::DistanceHistogram> model);

        protected:
            SimpleDataset data;
            std::unique_ptr<hist::DistanceHistogram> model;

            /**
             * @brief Reevaluate the model intensity, and prepare it for fitting. 
             */
            void refresh_model();

			/**
			 * @brief Splice values from the model to match the data.
			 * 
			 * @param ym the model y-values corresponding to xm
			 */
			[[nodiscard]] std::vector<double> splice(const std::vector<double>& ym) const;
    };
}