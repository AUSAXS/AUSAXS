#pragma once

#include <math/MathConcepts.h>
#include <math/Decomposition.h>
#include <math/Vector.h>
#include <math/MathFwd.h>

#include <memory>

namespace ausaxs {
    class LUPDecomposition : public Decomposition {
        public: 
            LUPDecomposition(const Matrix<double>& A);

            ~LUPDecomposition() override = default;

            // follows the C implementation from Wikipedia: https://en.wikipedia.org/wiki/LU_decomposition
            void decompose() override;

            double determinant() const;

            int permutations;
        private: 
            Vector<double> P;
            ::std::unique_ptr<Matrix<double>> Ap;
    };
}