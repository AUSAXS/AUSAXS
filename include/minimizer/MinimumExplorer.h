#pragma once

#include <minimizer/Minimizer.h>

namespace mini {
	/**
	 * @brief A class to explore the area near a minimum. 
     *        This is useful since most of our problems are not continuous and have varying step sizes. This class ensures that the function actually varies in the area explored.
	 */
	class MinimumExplorer : public Minimizer {
		public:
            MinimumExplorer(double(&func)(const double*), unsigned int evals);

            MinimumExplorer(std::function<double(const double*)> func, unsigned int evals);

            MinimumExplorer(double(&func)(const double*), const Parameter& param, unsigned int evals);

            MinimumExplorer(std::function<double(const double*)> func, const Parameter& param, unsigned int evals);

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
            Dataset2D landscape(unsigned int evals) override;

            /**
             * @brief Get the evaluated points and their function values.
             */
            Dataset2D get_evaluated_points() const override;

        private:
            unsigned int evals;

			/**
			 * @brief Perform the minimization.
			 */
			Result minimize_override() override;
	};
}