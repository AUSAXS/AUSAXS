#include "Matrix.h"
#include "Vector.h"
#include "LinearSolver.h"

#include <math.h>
#include <iostream>

class GivensSolver : public LinearSolver {
public:
    GivensSolver(const Matrix& A) : N(A.N), M(A.M) {decomp(A);}
    ~GivensSolver() override {}

    Vector solve(const Vector& b) const override {
        double theta, xp, xq;
        Vector x = b.copy();

        // start by applying the rotations to b
        for (int p = 0; p < M; p++) {
            for (int q = p+1; q < N; q++) {
                theta = G[q][p];
                xp = cos(theta)*x[p] + sin(theta)*x[q];
                xq = cos(theta)*x[q] - sin(theta)*x[p];
                x[p] = xp;
                x[q] = xq;
            }
        }

        // then make a back substitution to actually solve it
        for (int i = M-1; i >= 0; i--) {
            for (int j = i+1; j < M; j++)
                x[i] -= G[i][j]*x[j];
            x[i] /= G[i][i];
        }
        return x;
    }

private: 
    Matrix G;
    const int N, M;

    void decomp(const Matrix& A) {
        G = A.copy();
        double Aqi, Api;
        for (int p = 0; p < M; p++) {
            for (int q = p+1; q < N; q++) {
                double theta = atan2(A[q][p], A[p][p]);
                for (int i = p; i < M; i++) {
                    Api = A[p][i];
                    Aqi = A[q][i];
                    G[p][i] = Api*cos(theta) + Aqi*sin(theta);
                    G[q][i] = Aqi*cos(theta) - Api*sin(theta);
                }
                G[q][p] = theta; // replace with the rotation angle
            }
        }
    }
};