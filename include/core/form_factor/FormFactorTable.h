// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <constants/Constants.h>

#include <array>

// Five-Gaussian form factor table values. See each nested namespace for the source.
namespace ausaxs::constants::form_factor {
    namespace {
        constexpr double s_to_q_factor = 1./(4*4*std::numbers::pi*std::numbers::pi); // q = 4πs --> s = q/(4pi)

        /**
         * @brief Convert a Gaussian form factor from s to q.
         *        This is purely for convenience, such that the tabulated values are easier to read.
         */
        constexpr std::array<double, 5> s_to_q(std::array<double, 5> a) {
            for (unsigned int i = 0; i < 5; ++i) {
                a[i] *= s_to_q_factor;
            }
            return a;
        }
    }

    // International Tables for Crystallography, https://lampx.tugraz.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php        
    namespace H {
        constexpr std::array<double, 5> a = {0.489918, 0.262003, 0.196767, 0.049879, 0};
        constexpr std::array<double, 5> b = s_to_q({20.6593, 7.74039, 49.5519, 2.20159, 0});
        constexpr double c = 0.001305;
    }

    // Waasmeier & Kirfel, https://doi.org/10.1107/S0108767394013292
    namespace O {
        constexpr std::array<double, 5> a =        { 2.960427, 2.508818, 0.637853,  0.722838, 1.142756};
        constexpr std::array<double, 5> b = s_to_q({14.182259, 5.936858, 0.112726, 34.958481, 0.390240});
        constexpr double c = 0.027014;
    }

    // Waasmeier & Kirfel, https://doi.org/10.1107/S0108767394013292
    namespace N {
        constexpr std::array<double, 5> a =        {11.893780,  3.277479,  1.858092, 0.858927, 0.912985};
        constexpr std::array<double, 5> b = s_to_q({ 0.000158, 10.232723, 30.344690, 0.656065, 0.217287});
        constexpr double c = -11.804902;
    }

    // Waasmeier & Kirfel, https://doi.org/10.1107/S0108767394013292
    namespace S {
        constexpr std::array<double, 5> a =        {6.362157, 5.154568,  1.473732, 1.635073,  1.209372};
        constexpr std::array<double, 5> b = s_to_q({1.514347, 22.092528, 0.061373, 55.445176, 0.646925});
        constexpr double c = 0.154722;
    }

    // Waasmeier & Kirfel, https://doi.org/10.1107/S0108767394013292
    namespace C {
        constexpr std::array<double, 5> a =        { 2.657506, 1.078079,  1.490909, -4.241070, 0.713791};
        constexpr std::array<double, 5> b = s_to_q({14.780758, 0.776775, 42.086843, -0.000294, 0.239535});
        constexpr double c = 4.297983;
    }

    // Waasmeier & Kirfel, https://doi.org/10.1107/S0108767394013292
    // this is just the form factor of argon
    namespace other {
        constexpr std::array<double, 5> a =        {7.188004, 6.638454,  0.454180,  1.929593,  1.523654};
        constexpr std::array<double, 5> b = s_to_q({0.956221, 15.339877, 15.339862, 39.043824, 0.062409});
        constexpr double c = 0.265954;
    }

    namespace excluded_volume {
        constexpr std::array<double, 5> a = {1, 0, 0, 0, 0};
        constexpr std::array<double, 5> b = {radius::average_atomic_radius*radius::average_atomic_radius/2, 0, 0, 0, 0};
        constexpr double c = 0;
    }

    // Grudinin, Garkavenko, & Kazennov, https://doi.org/10.1107/s2059798317005745
    namespace CH_sp3 {
        constexpr std::array<double, 5> a =        { 2.909530, 0.485267,  1.516151,  0.206905, 1.541626};
        constexpr std::array<double, 5> b = s_to_q({13.933084, 23.221524, 41.990403, 4.974183, 0.679266});
        constexpr double c = 0.337670;
    }

    // Grudinin, Garkavenko, & Kazennov, https://doi.org/10.1107/s2059798317005745
    namespace CH2_sp3 {
        constexpr std::array<double, 5> a =        { 3.275723, 0.870037,  1.534606,  0.395078, 1.544562};
        constexpr std::array<double, 5> b = s_to_q({13.408502, 23.785175, 41.922444, 5.019072, 0.724439});
        constexpr double c = 0.377096;
    }

    // Grudinin, Garkavenko, & Kazennov, https://doi.org/10.1107/s2059798317005745
    namespace CH3_sp3 {
        constexpr std::array<double, 5> a =        { 3.681341, 1.228691,  1.549320,  0.574033, 1.554377};
        constexpr std::array<double, 5> b = s_to_q({13.026207, 24.131974, 41.869426, 4.984373, 0.765769});
        constexpr double c = 0.409294;
    }

    // Grudinin, Garkavenko, & Kazennov, https://doi.org/10.1107/s2059798317005745
    namespace CH_sp2 {
        constexpr std::array<double, 5> a =        { 2.909457, 0.484873,  1.515916,  0.207091, 1.541518};
        constexpr std::array<double, 5> b = s_to_q({13.934162, 23.229153, 41.991425, 4.983276, 0.679898});
        constexpr double c = 0.338296;
    }

    // Grudinin, Garkavenko, & Kazennov, https://doi.org/10.1107/s2059798317005745
    namespace CH_arom {
        constexpr std::array<double, 5> a =        { 2.168070, 1.275811,  1.561096,  0.742395, -6.151144};
        constexpr std::array<double, 5> b = s_to_q({12.642907, 18.420069, 41.768517, 1.535360, -0.045937});
        constexpr double c = 7.400917;
    }

    // Grudinin, Garkavenko, & Kazennov, https://doi.org/10.1107/s2059798317005745
    namespace OH_alc {
        constexpr std::array<double, 5> a =        { 0.456221, 3.219608,  0.812773,  2.666928, 1.380927};
        constexpr std::array<double, 5> b = s_to_q({21.503498, 13.397134, 34.547137, 5.826620, 0.412902});
        constexpr double c = 0.463202;
    }

    // Grudinin, Garkavenko, & Kazennov, https://doi.org/10.1107/s2059798317005745
    namespace OH_acid {
        constexpr std::array<double, 5> a =        { 3.213280, 0.463019,  0.815724,  2.664450, 1.384266};
        constexpr std::array<double, 5> b = s_to_q({13.383078, 21.362223, 34.531415, 5.823549, 0.410805});
        constexpr double c = 0.458919;
    }

    // Grudinin, Garkavenko, & Kazennov, https://doi.org/10.1107/s2059798317005745
    namespace O_res {
        constexpr std::array<double, 5> a =        { 0.688944, 2.929687, 0.416472,  2.606983,  1.319232};
        constexpr std::array<double, 5> b = s_to_q({29.319200, 6.572228, 64.951658, 16.267799, 0.455640});
        constexpr double c = 0.537548;
    }

    // Grudinin, Garkavenko, & Kazennov, https://doi.org/10.1107/s2059798317005745
    namespace NH {
        constexpr std::array<double, 5> a =        { 1.650531, 0.429639, 2.144736,  1.851894,  1.408921};
        constexpr std::array<double, 5> b = s_to_q({10.603730, 6.987283, 29.939901, 10.573859, 0.611678});
        constexpr double c = 0.510589;
    }

    // Grudinin, Garkavenko, & Kazennov, https://doi.org/10.1107/s2059798317005745
    namespace NH2 {
        constexpr std::array<double, 5> a =        { 1.904157, 1.942536,  2.435585,  0.730512, 1.379728};
        constexpr std::array<double, 5> b = s_to_q({10.803702, 10.792421, 29.610479, 6.847755, 0.709687});
        constexpr double c = 0.603738;
    }

    // Grudinin, Garkavenko, & Kazennov, https://doi.org/10.1107/s2059798317005745
    namespace NH_plus {
        constexpr std::array<double, 5> a =        { 1.426540, 0.426903, 1.878894,  1.608251,  1.200216};
        constexpr std::array<double, 5> b = s_to_q({10.652268, 7.017651, 29.878525, 10.619493, 0.631765});
        constexpr double c = 0.456024;
    }

    // Grudinin, Garkavenko, & Kazennov, https://doi.org/10.1107/s2059798317005745
    namespace NH2_plus {
        constexpr std::array<double, 5> a =        { 3.823896, 0.531490,  1.713620,  0.322552, 1.287502};
        constexpr std::array<double, 5> b = s_to_q({10.305028, 25.631593, 30.215026, 3.576178, 0.506824});
        constexpr double c = 0.317728;
    }

    // Grudinin, Garkavenko, & Kazennov, https://doi.org/10.1107/s2059798317005745
    namespace NH3_plus {
        constexpr std::array<double, 5> a =        { 1.882162, 1.933200,  2.465843,  0.927311, 1.190889};
        constexpr std::array<double, 5> b = s_to_q({10.975157, 10.956008, 29.208572, 6.663555, 0.843650});
        constexpr double c = 0.597322;
    }

    // Grudinin, Garkavenko, & Kazennov, https://doi.org/10.1107/s2059798317005745
    namespace NH_guanine {
        constexpr std::array<double, 5> a =        { 3.630164, 0.228310,  1.869734,  0.170550, 1.440894};
        constexpr std::array<double, 5> b = s_to_q({10.267139, 25.118086, 30.241288, 3.412776, 0.486644});
        constexpr double c = 0.323504;
    }

    // Grudinin, Garkavenko, & Kazennov, https://doi.org/10.1107/s2059798317005745
    namespace NH2_guanine {
        constexpr std::array<double, 5> a =        { 1.792216, 0.724464, 2.347044,  1.903020,  1.313042};
        constexpr std::array<double, 5> b = s_to_q({10.830060, 6.846763, 29.579607, 10.800018, 0.720448});
        constexpr double c = 0.583312;
    }

    // Grudinin, Garkavenko, & Kazennov, https://doi.org/10.1107/s2059798317005745
    namespace SH {
        constexpr std::array<double, 5> a =        { 0.570042, 6.337416, 1.641643,  5.398549,  1.527982};
        constexpr std::array<double, 5> b = s_to_q({11.447986, 1.197657, 55.401032, 22.420955, 2.356552});
        constexpr double c = 1.523944;
    }
}