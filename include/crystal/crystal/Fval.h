// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <crystal/miller/Miller.h>
#include <math/Vector3.h>
#include <data/DataFwd.h>
#include <utility/UtilityFwd.h>

#include <vector>
#include <complex>

namespace ausaxs::crystal {
    class Fval {
        public: 
            Fval();
            Fval(int h, int k, int l);

            static void set_points(std::vector<Vector3<double>>&& points);

            static void set_basis(const Basis3D& basis);

            static std::vector<Vector3<double>>& get_points();

            static Basis3D get_basis();

            static Vector3<double> Q(const Miller& miller);

            double I() const;

            static data::Molecule as_protein();
            
            Miller hkl;
            double qlength;
            Vector3<double> q;
            std::complex<double> fval;
        private: 
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

            static void precompute_factors();
    };
}