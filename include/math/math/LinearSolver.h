#pragma once

#include <math/Matrix.h>
#include <math/Vector.h>

namespace ausaxs {
	class LinearSolver {
		public:
			virtual ~LinearSolver(){}

			enum Algorithm {Givens};

			/**
			 * @brief Solve a linear equation of the form Ax = b through a QR decomposition.
			 */
			virtual Vector<double> solve(const Vector<double>& b) const = 0;

		private: 
			Algorithm algorithm;
	};
}