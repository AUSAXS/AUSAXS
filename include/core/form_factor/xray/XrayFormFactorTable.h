// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <constants/Constants.h>

#include <array>
#include <numbers>

// Five-Gaussian X-ray form factor table values. See each entry for its source.
namespace ausaxs::form_factor::xray::coefficients {
    /**
     * @brief The coefficients of a five-Gaussian form factor approximation, f(q) = sum_i a_i exp(-b_i q^2) + c.
     *        The b coefficients are stored in q-units (Å^2).
     */
    struct FiveGaussian {
        std::array<double, 5> a;
        std::array<double, 5> b;
        double c;
    };

    namespace {
        constexpr double s_to_q_factor = 1./(4*4*std::numbers::pi*std::numbers::pi); // q = 4πs --> s = q/(4pi)

        /**
         * @brief Convert a Gaussian form factor from s to q.
         *        This is purely for convenience, such that the tabulated values are easier to read.
         */
        constexpr std::array<double, 5> s_to_q(std::array<double, 5> a) {
            for (int i = 0; i < 5; ++i) {
                a[i] *= s_to_q_factor;
            }
            return a;
        }
    }

    // International Tables for Crystallography, https://lampx.tugraz.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php
    constexpr FiveGaussian H           {.a = { 0.489918, 0.262003, 0.196767,  0.049879,         0}, .b = s_to_q({  20.6593,    7.74039,    49.5519,    2.20159,         0}), .c =   0.001305};

    // Waasmeier & Kirfel, https://doi.org/10.1107/S0108767394013292
    constexpr FiveGaussian C           {.a = { 2.657506, 1.078079, 1.490909, -4.241070,  0.713791}, .b = s_to_q({14.780758,   0.776775,  42.086843,  -0.000294,  0.239535}), .c =   4.297983};
    constexpr FiveGaussian N           {.a = {11.893780, 3.277479, 1.858092,  0.858927,  0.912985}, .b = s_to_q({ 0.000158,  10.232723,  30.344690,   0.656065,  0.217287}), .c = -11.804902};
    constexpr FiveGaussian O           {.a = { 2.960427, 2.508818, 0.637853,  0.722838,  1.142756}, .b = s_to_q({14.182259,   5.936858,   0.112726,  34.958481,  0.390240}), .c =   0.027014};
    constexpr FiveGaussian F           {.a = { 3.511943, 2.772244, 0.678385,  0.915159,  1.089261}, .b = s_to_q({10.687859,   4.380466,   0.093982,  27.255203,  0.313066}), .c =   0.032557};
    constexpr FiveGaussian Na          {.a = { 4.910127, 3.081783, 1.262067,  1.098938,  0.560991}, .b = s_to_q({ 3.281434,   9.119178,   0.102763, 132.013942,  0.405878}), .c =   0.079712};
    constexpr FiveGaussian Mg          {.a = { 4.708971, 1.194814, 1.558157,  1.170413,  3.239403}, .b = s_to_q({ 4.875207, 108.506079,   0.111516,  48.292407,  1.928171}), .c =   0.126842};
    constexpr FiveGaussian P           {.a = { 1.950541, 4.146930, 1.494560,  1.522042,  5.729711}, .b = s_to_q({ 0.908139,  27.044953,   0.071280,  67.520190,  1.981173}), .c =   0.155233};
    constexpr FiveGaussian S           {.a = { 6.372157, 5.154568, 1.473732,  1.635073,  1.209372}, .b = s_to_q({ 1.514347,  22.092528,   0.061373,  55.445176,  0.646925}), .c =   0.154722};
    constexpr FiveGaussian Cl          {.a = { 1.446071, 6.870609, 6.151801,  1.750347,  0.634168}, .b = s_to_q({ 0.052357,   1.193165,  18.343416,  46.398394,  0.401005}), .c =   0.146773};
    constexpr FiveGaussian Ar          {.a = { 7.188004, 6.638454, 0.454180,  1.929593,  1.523654}, .b = s_to_q({ 0.956221,  15.339877,  15.339862,  39.043824,  0.062409}), .c =   0.265954}; 
    constexpr FiveGaussian K           {.a = { 8.163991, 7.146945, 1.070140,  0.877316,  1.486434}, .b = s_to_q({12.816323,   0.808945, 210.327009,  39.597651,  0.052821}), .c =   0.253614};
    constexpr FiveGaussian Ca          {.a = { 8.593655, 1.477324, 1.436254,  1.182839,  7.113258}, .b = s_to_q({10.460644,   0.041891,  81.390382, 169.847839,  0.688098}), .c =   0.196255};
    constexpr FiveGaussian Mn          {.a = {11.709542, 1.733414, 2.673141,  2.023368,  7.003180}, .b = s_to_q({ 5.597120,   0.017800,  21.788419,  89.517915,  0.383054}), .c =  -0.147293};
    constexpr FiveGaussian Fe          {.a = {12.311098, 1.876623, 3.066177,  2.070451,  6.975185}, .b = s_to_q({ 5.009415,   0.014461,  18.743041,  82.767874,  0.346506}), .c =  -0.304931};
    constexpr FiveGaussian Co          {.a = {12.914510, 2.481908, 3.466894,  2.106351,  6.960892}, .b = s_to_q({ 4.507138,   0.009126,  16.438130,  76.987317,  0.314418}), .c =  -0.936572};
    constexpr FiveGaussian Ni          {.a = {13.521865, 6.947285, 3.866028,  2.135900,  4.284731}, .b = s_to_q({ 4.077277,   0.286763,  14.622634,  71.966078,  0.004437}), .c =  -2.762697};
    constexpr FiveGaussian Cu          {.a = {14.014192, 4.784577, 5.056806,  1.457971,  6.932996}, .b = s_to_q({ 3.738280,   0.003744,  13.034982,  72.554793,  0.265666}), .c =  -3.254477};
    constexpr FiveGaussian Zn          {.a = {14.741002, 6.907748, 4.642337,  2.191766, 38.424042}, .b = s_to_q({ 3.388232,   0.243315,  11.903689,  63.312130,  0.000397}), .c = -36.915828};
    constexpr FiveGaussian Se          {.a = {17.354071, 4.653248, 4.259489,  4.136455,  6.749163}, .b = s_to_q({ 2.349787,   0.002550,  15.579460,  45.181201,  0.177432}), .c =  -3.160982};
    constexpr FiveGaussian Br          {.a = {17.550570, 5.411882, 3.937180,  3.880645,  6.707793}, .b = s_to_q({ 2.119226,  16.557185,   0.002481,  42.164009,  0.162121}), .c =  -2.492088};
    constexpr FiveGaussian I           {.a = {19.884502, 6.736593, 8.110516,  1.170953, 17.548715}, .b = s_to_q({ 4.628591,   0.027754,  31.849096,  84.406391,  0.463550}), .c =  -0.448811};
    constexpr FiveGaussian other = Ar; // default form factor for unknown atoms

    // Grudinin, Garkavenko, & Kazennov, https://doi.org/10.1107/s2059798317005745
    constexpr FiveGaussian CH_sp3      {.a = { 2.909530, 0.485267, 1.516151,  0.206905,  1.541626}, .b = s_to_q({13.933084,  23.221524,  41.990403,   4.974183,  0.679266}), .c =   0.337670};
    constexpr FiveGaussian CH2_sp3     {.a = { 3.275723, 0.870037, 1.534606,  0.395078,  1.544562}, .b = s_to_q({13.408502,  23.785175,  41.922444,   5.019072,  0.724439}), .c =   0.377096};
    constexpr FiveGaussian CH3_sp3     {.a = { 3.681341, 1.228691, 1.549320,  0.574033,  1.554377}, .b = s_to_q({13.026207,  24.131974,  41.869426,   4.984373,  0.765769}), .c =   0.409294};
    constexpr FiveGaussian CH_sp2      {.a = { 2.909457, 0.484873, 1.515916,  0.207091,  1.541518}, .b = s_to_q({13.934162,  23.229153,  41.991425,   4.983276,  0.679898}), .c =   0.338296};
    constexpr FiveGaussian CH_arom     {.a = { 2.168070, 1.275811, 1.561096,  0.742395, -6.151144}, .b = s_to_q({12.642907,  18.420069,  41.768517,   1.535360, -0.045937}), .c =   7.400917};
    constexpr FiveGaussian OH_alc      {.a = { 0.456221, 3.219608, 0.812773,  2.666928,  1.380927}, .b = s_to_q({21.503498,  13.397134,  34.547137,   5.826620,  0.412902}), .c =   0.463202};
    constexpr FiveGaussian OH_acid     {.a = { 3.213280, 0.463019, 0.815724,  2.664450,  1.384266}, .b = s_to_q({13.383078,  21.362223,  34.531415,   5.823549,  0.410805}), .c =   0.458919};
    constexpr FiveGaussian O_res       {.a = { 0.688944, 2.929687, 0.416472,  2.606983,  1.319232}, .b = s_to_q({29.319200,   6.572228,  64.951658,  16.267799,  0.455640}), .c =   0.537548};
    constexpr FiveGaussian NH          {.a = { 1.650531, 0.429639, 2.144736,  1.851894,  1.408921}, .b = s_to_q({10.603730,   6.987283,  29.939901,  10.573859,  0.611678}), .c =   0.510589};
    constexpr FiveGaussian NH2         {.a = { 1.904157, 1.942536, 2.435585,  0.730512,  1.379728}, .b = s_to_q({10.803702,  10.792421,  29.610479,   6.847755,  0.709687}), .c =   0.603738};
    constexpr FiveGaussian NH_plus     {.a = { 1.426540, 0.426903, 1.878894,  1.608251,  1.200216}, .b = s_to_q({10.652268,   7.017651,  29.878525,  10.619493,  0.631765}), .c =   0.456024};
    constexpr FiveGaussian NH2_plus    {.a = { 3.823896, 0.531490, 1.713620,  0.322552,  1.287502}, .b = s_to_q({10.305028,  25.631593,  30.215026,   3.576178,  0.506824}), .c =   0.317728}; // NOLINT(modernize-use-std-numbers): tabulated value, not 1/pi
    constexpr FiveGaussian NH3_plus    {.a = { 1.882162, 1.933200, 2.465843,  0.927311,  1.190889}, .b = s_to_q({10.975157,  10.956008,  29.208572,   6.663555,  0.843650}), .c =   0.597322};
    constexpr FiveGaussian NH_guanine  {.a = { 3.630164, 0.228310, 1.869734,  0.170550,  1.440894}, .b = s_to_q({10.267139,  25.118086,  30.241288,   3.412776,  0.486644}), .c =   0.323504};
    constexpr FiveGaussian NH2_guanine {.a = { 1.792216, 0.724464, 2.347044,  1.903020,  1.313042}, .b = s_to_q({10.830060,   6.846763,  29.579607,  10.800018,  0.720448}), .c =   0.583312};
    constexpr FiveGaussian SH          {.a = { 0.570042, 6.337416, 1.641643,  5.398549,  1.527982}, .b = s_to_q({11.447986,   1.197657,  55.401032,  22.420955,  2.356552}), .c =   1.523944};

    constexpr FiveGaussian excluded_volume {.a = {1, 0, 0, 0, 0}, .b = {constants::radius::average_atomic_radius*constants::radius::average_atomic_radius/2, 0, 0, 0, 0}, .c = 0};
}
