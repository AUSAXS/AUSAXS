// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <mini/Minimizer.h>

namespace ausaxs::mini {
	/**
	 * @brief A scanning minimizer. 
	 * 		  This algorithm performs a simple scan of the chi2 landscape, and then uses a local minimizer to find the minimum. 
	 * 		  This should only be used for relatively simple optimization problems. 
	 */
	class Scan : public Minimizer {
		public:
			Scan() = default;

			Scan(double(&func)(std::vector<double>), unsigned int evals = 100);

            Scan(std::function<double(std::vector<double>)> func, unsigned int evals = 100);

            Scan(double(&func)(std::vector<double>), const Parameter& param, unsigned int evals = 100);

            Scan(std::function<double(std::vector<double>)> func, const Parameter& param, unsigned int evals = 100);

            /**
             * @brief Destructor.
             */
            virtual ~Scan() override = default;

            /**
             * @brief Generate a landscape of the function.
             */
            mini::Landscape landscape(unsigned int evals) override;

			void add_parameter(const Parameter& param) override;

		private:
			/**
			 * @brief Perform the minimization.
			 */
			Result minimize_override() override;
	};
}