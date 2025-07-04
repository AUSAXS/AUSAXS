// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/Matrix.h>
#include <math/Vector.h>
#include <math/LinearSolver.h>

#include <math.h>

namespace ausaxs {
    class GivensSolver : public LinearSolver {
        public:
            // GivensSolver(const Matrix<double>& A) : N(A.N), M(A.M) {decomp(A);}
            GivensSolver(const Matrix<double>& A) = delete;
            ~GivensSolver() override = default;

            Vector<double> solve(const Vector<double>& b) const override {
                double theta, xp, xq;
                Vector<double> x = b;

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
            Matrix<double> G;
            const int N, M;

            void decomp(const Matrix<double>& A) {
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
}