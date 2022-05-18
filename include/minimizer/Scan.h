#pragma once

#include <minimizer/Minimizer.h>

namespace mini {
	/**
	 * @brief A scanning minimizer. This algorithm performs scans of ever-increasing resolution until a minimum is found. This should only be used for nice and easy optimization problems. 
	 */
	class Scan : public Minimizer {
		public:
			/**
			 * @brief Perform the minimization.
			 */
			Result minimize() override;

			/**
			 * @brief Add a parameter.
			 * 
			 * @param par The name of the parameter.
			 * @param guess The start value of the parameter. 
			 */
			void add_parameter(std::string par, double guess);

			/**
			 * @brief Add a parameter with limits.
			 * 
			 * @param par The name of the parameter.
			 * @param guess The start value of the parameter. 
			 * @param limits The limits of the parameter.
			 */
			void add_parameter(std::string par, double guess, Limit limits);
	};
}