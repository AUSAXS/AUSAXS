#include "Decomposition.h"
#include "Matrix.h"
#include "Vector.h"
#include "Slice.h"

// QR decomposition by Gram-Schmidt orthogonalization
class QRDecomposition : public Decomposition {
public: 
    QRDecomposition(const Matrix& A) : Q(A) {decompose();}

    Matrix inverse() const {
		// basically we just solve m equations of the form Ax = e_i, and construct A^-1 from the m solutions to this equation
		Matrix A(Q.M, Q.N);
		Vector e(Q.N);
		for (size_t i = 0; i < Q.N; i++) {
			e[i] = 1;
			A.col(i) = solve(e);
			e[i] = 0;
		}
		return A;
	}

	// solve an equation of the form QRx = b
	Vector solve(Vector b) const {
		b = Q.T()*b;
		Vector x(b.N);
		for (int i = x.N-1; i >= 0; i--) {
			x[i] = b[i];
			for (size_t j = i+1; j < b.N; j++)
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
        R = Matrix(Q.N, Q.M);
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

    Matrix Q, R;
};