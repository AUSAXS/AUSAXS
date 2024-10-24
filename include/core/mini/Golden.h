#pragma once

#include <mini/Minimizer.h>

namespace ausaxs::mini {
	/**
	 * @brief A bisection algorithm based on the golden ratio. This only supports unitary problems. 
	 */
	class Golden : public Minimizer {
		public:
            Golden() = default;

            Golden(double(&func)(std::vector<double>));

            Golden(std::function<double(std::vector<double>)> func);

            Golden(double(&func)(std::vector<double>), const Parameter& param);

            Golden(std::function<double(std::vector<double>)> func, const Parameter& param);

            /**
             * @brief Destructor.
             */
            ~Golden() override = default;

			/**
			 * @brief Add a parameter.
			 */
			void add_parameter(const Parameter& param) override;

        private:
            inline static const double phi = (1 + std::sqrt(5))/2;
            inline static const double invphi = 1/phi;
            inline static const double invphi2 = invphi*invphi;

			/**
			 * @brief Perform the minimization.
			 */
			Result minimize_override() override;

            /**
             * @brief Golden-section search. 
             * 
             * Given a function with a single local minimum in the interval [a, b], this search finds and returns the interval [c, d] where d-c < tol.
             * 
			 * @param bounds The bounds to search within. 
             */
            Limit search(Limit bounds) const;
	};
}