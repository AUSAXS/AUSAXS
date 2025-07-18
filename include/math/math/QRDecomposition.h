// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/Decomposition.h>
#include <math/Matrix.h>
#include <math/Vector.h>
#include <math/slices/Slice.h>

namespace ausaxs {
	// QR decomposition by Gram-Schmidt orthogonalization
	class QRDecomposition : public Decomposition {
		public: 
			QRDecomposition(const Matrix<double>& A) : Q(A) {decompose();}

			Matrix<double> inverse() const {
				// basically we just solve m equations of the form Ax = e_i, and construct A^-1 from the m solutions to this equation
				Matrix<double> A(Q.M, Q.N);
				Vector<double> e(Q.N);
				for (size_t i = 0; i < Q.N; i++) {
					e[i] = 1;
					A.col(i) = solve(e);
					e[i] = 0;
				}
				return A;
			}

			// solve an equation of the form QRx = b
			Vector<double> solve(Vector<double> b) const {
				b = Q.T()*b;
				Vector<double> x(b.size());
				for (int i = x.size()-1; i >= 0; i--) {
					x[i] = b[i];
					for (size_t j = i+1; j < b.size(); j++)
						x[i] -= R[i][j]*x[j];
					x[i] /= R[i][i];
				}
				return x;
			}

			// up to a sign
			double abs_determinant() const {
				double det = R[0][0];
				for (size_t i = 1; i < R.N; i++) {
					det *= R[i][i];
				}
				return det;
			}

			void decompose() override {
				R = Matrix<double>(Q.N, Q.M);
				double ujvi, ujuj;
				for (size_t i = 0; i < Q.M; i++) {
					for (size_t j = 0; j < i; j++) {
						ujvi = Q.col(j).dot(Q.col(i)); 
						ujuj = Q.col(j).dot(Q.col(j));
						R[j][i] = ujvi/ujuj; // a_ij
						Q.col(i) -= Q.col(j)*R[j][i];
						R[j][i] *= R[j][j];
					}
					R[i][i] = Q.col(i).norm();
				}
				for (size_t i = 0; i < Q.M; i++) {
					Q.col(i) = Q.col(i)/R[i][i];
				}
			}

			Matrix<double> Q, R;
	};
}