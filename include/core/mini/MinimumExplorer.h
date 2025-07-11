// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <mini/Minimizer.h>

namespace ausaxs::mini {
	/**
	 * @brief A class to explore the area near a minimum. 
     *        This is useful since most of our problems are not continuous and have varying step sizes. This class ensures that the function actually varies in the area explored.
	 */
	class MinimumExplorer : public Minimizer {
		public:
            MinimumExplorer() = default;

            MinimumExplorer(double(&func)(std::vector<double>), unsigned int evals = 100);

            MinimumExplorer(std::function<double(std::vector<double>)> func, unsigned int evals = 100);

            MinimumExplorer(double(&func)(std::vector<double>), const Parameter& param, unsigned int evals = 100);

            MinimumExplorer(std::function<double(std::vector<double>)> func, const Parameter& param, unsigned int evals = 100);

            /**
             * @brief Destructor.
             */
            ~MinimumExplorer() override = default;

			/**
			 * @brief Add a parameter.
			 */
			void add_parameter(const Parameter& param) override;

            /**
             * @brief Generate a landscape of the function.
             */
            mini::Landscape landscape(unsigned int evals = 100) override;

        private:
			/**
			 * @brief Perform the minimization.
			 */
			Result minimize_override() override;
	};
}