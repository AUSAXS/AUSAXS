#pragma once

#include <crystal/miller/Miller.h>
#include <math/Vector3.h>

#include <vector>
#include <complex>

namespace crystal {
    class Fval {
        public: 
            Fval();
            Fval(double h, double k, double l);

            static void set_points(std::vector<Vector3<double>>&& points);

            static std::vector<Vector3<double>>& get_points();

            double I() const;

            Miller<> hkl;
            double qlength;
            Vector3<double> q;
            std::complex<double> fval;
            Vector3<double> ap, bp, cp;
        private: 
            std::complex<double> i = std::complex<double>(0, 1);
            inline static std::vector<Vector3<double>> points;

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