#pragma once

#include <minimizer/Minimizer.h>

namespace mini {
	/**
	 * @brief A bisection algorithm based on the golden ratio. This only supports unitary problems. 
	 */
	class Golden : public Minimizer {
		public:
            Golden(double(&func)(double*), std::string par, Limit bounds);

            Golden(std::function<double(double*)> func, std::string par, Limit bounds);

			/**
			 * @brief Perform the minimization.
			 */
			Result minimize() const override;

			/**
			 * @brief Add a parameter.
			 * 
			 * @param par The name of the parameter.
			 * @param guess The start value of the parameter. 
			 */
			void add_parameter(std::string par, double guess) override;

			/**
			 * @brief Add a parameter with limits.
			 * 
			 * @param par The name of the parameter.
			 * @param guess The start value of the parameter. 
			 * @param bounds The bounds of the parameter.
			 */
			void add_parameter(std::string par, double guess, Limit bounds) override;

			/**
			 * @brief Add a parameter with limits.
			 * 
			 * @param par The name of the parameter.
			 * @param bounds The bounds on the parameter.
			 */
			void add_parameter(std::string par, Limit bounds);

            /**
             * @brief Generate a landscape of the function.
             */
            Dataset landscape(unsigned int evals = 100) const;

            /**
             * @brief Get the evaluated points and their function values.
             */
            Dataset get_evaluated_points() const;

        private:
            inline static const double phi = (1 + std::sqrt(5))/2;
            inline static const double invphi = 1/phi;
            inline static const double invphi2 = invphi*invphi;
            Limit bounds;

            /**
             * @brief Golden-section search. 
             * 
             * Given a function with a single local minimum in the interval [a, b], this search finds and returns the interval [c, d] where d-c < tol.
             * 
             * @param a Lower bound containing the solution.
             * @param b Upper bound containing the solution.
             */
            Limit search(double a, double b) const;
	};
}