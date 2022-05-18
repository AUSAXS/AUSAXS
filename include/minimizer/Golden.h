#pragma once

#include <minimizer/Minimizer.h>

namespace mini {
	/**
	 * @brief A bisection algorithm based on the golden ratio. This only supports unitary problems. 
	 */
	class Golden : public Minimizer {
		public:
            Golden(double(&func)(double*), const Parameter& param);

            Golden(std::function<double(double*)> func, const Parameter& param);

			/**
			 * @brief Perform the minimization.
			 */
			Result minimize() const override;

			/**
			 * @brief Add a parameter.
			 */
			void add_parameter(const Parameter& param) override;

            /**
             * @brief Generate a landscape of the function.
             */
            Dataset landscape(unsigned int evals = 100) const override;

            /**
             * @brief Get the evaluated points and their function values.
             */
            Dataset get_evaluated_points() const;

        private:
            inline static const double phi = (1 + std::sqrt(5))/2;
            inline static const double invphi = 1/phi;
            inline static const double invphi2 = invphi*invphi;

            /**
             * @brief Golden-section search. 
             * 
             * Given a function with a single local minimum in the interval [a, b], this search finds and returns the interval [c, d] where d-c < tol.
             * 
			 * @param bounds The bounds to search within. 
             */
            Limit search(const Limit bounds) const;
	};
}