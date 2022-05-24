#pragma once

#include <minimizer/Minimizer.h>

namespace mini {
	/**
	 * @brief A scanning minimizer. This algorithm performs scans of ever-increasing resolution until a minimum is found. This should only be used for nice and easy optimization problems. 
	 */
	class Scan : public Minimizer {
		public:
            /**
             * @brief Generate a landscape of the function.
             */
            Dataset landscape(unsigned int evals) const override;

            /**
             * @brief Get the evaluated points and their function values.
             */
            Dataset get_evaluated_points() const override;

			void add_parameter(const Parameter& param);

		private:
			unsigned int bins = 100;

			/**
			 * @brief Perform the minimization.
			 */
			Result minimize_override() override;

			void looper(std::vector<double>& p, unsigned int index) const;
	};
}