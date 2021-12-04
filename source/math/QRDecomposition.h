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
			A[i] = solve(e);
			e[i] = 0;
		}
		return A.T();
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

	void decompose() override {
        R = Matrix(Q.N, Q.M);
		double ujvi, ujuj;
		for (size_t i = 0; i < Q.M; i++) {
			for (size_t j = 0; j < i; j++) {
				ujvi = Q[j].dot(Q[i]); 
                ujuj = Q[j].dot(Q[j]);
				R[i][j] = ujvi/ujuj; // a_ij
				Q[i] -= Q[j]*R[j][i];
				R[j][i] *= R[j][j];
			}
			R[i][i] = Q[i].norm();
		}
		for (size_t i = 0; i < Q.M; i++)
			Q[i] = Q[i]/R[i][i];
	}

	// void decompose() override {
    //     R = Matrix(Q.N, Q.M);
	// 	double ujvi, ujuj;
	// 	for (size_t i = 0; i < Q.M; i++) {
	// 		for (size_t j = 0; j < i; j++) {
	// 			ujvi = Q[j].dot(Q[i]); 
    //             ujuj = Q[j].dot(Q[j]);
	// 			R[j][i] = ujvi/ujuj; // a_ij
	// 			Q[i] -= Q[j]*R[j][i];
	// 			R[j][i] *= R[j][j];
	// 		}
	// 		R[i][i] = Q[i].norm();
	// 	}
	// 	for (size_t i = 0; i < Q.M; i++)
	// 		Q[i] = Q[i]/R[i][i];
	// }

    Matrix Q, R;
};