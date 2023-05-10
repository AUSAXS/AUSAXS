#pragma once

#include <crystal/miller/Miller.h>
#include <math/Vector3.h>

#include <vector>
#include <complex>

class Basis3D;
namespace crystal {
    class Fval {
        public: 
            Fval();
            Fval(double h, double k, double l);

            static void set_points(std::vector<Vector3<double>>&& points);

            static void set_basis(const Basis3D& basis);

            static std::vector<Vector3<double>>& get_points();

            static Basis3D get_basis();

            double I() const;

            Miller hkl;
            double qlength;
            Vector3<double> q;
            std::complex<double> fval;
        private: 
            std::complex<double> i = std::complex<double>(0, 1);
            inline static std::vector<Vector3<double>> points;
            inline static Vector3<double> ap, bp, cp;

            /**
             * @brief Calculate the F value for the given Miller indices
             */
            std::complex<double> F();

            /**
             * @brief Calculate the Q vector for the given Miller indices
             */
            Vector3<double> Q() const;
    };
}