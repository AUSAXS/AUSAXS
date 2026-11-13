#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/FormFactorInfo.h>
#include <form_factor/xray/XrayFormFactor.h>
#include <form_factor/xray/XrayFormFactorTable.h>

using namespace ausaxs;
using namespace form_factor;
using namespace form_factor::xray::coefficients;

// The five-Gaussian fits to the atomic groups which were used before the groups were built from their heavy atom and hydrogens.
// Grudinin, Garkavenko, & Kazennov, https://doi.org/10.1107/s2059798317005745
// These published values do not reproduce the curves Pepsi-SAXS itself computes from them, so they are only a rough reference.
namespace grudinin {
    constexpr FiveGaussian CH_sp3      {.a = { 2.909530, 0.485267, 1.516151,  0.206905,  1.541626}, .b = s_to_q({13.933084,  23.221524,  41.990403,   4.974183,  0.679266}), .c =   0.337670};
    constexpr FiveGaussian CH2_sp3     {.a = { 3.275723, 0.870037, 1.534606,  0.395078,  1.544562}, .b = s_to_q({13.408502,  23.785175,  41.922444,   5.019072,  0.724439}), .c =   0.377096};
    constexpr FiveGaussian CH3_sp3     {.a = { 3.681341, 1.228691, 1.549320,  0.574033,  1.554377}, .b = s_to_q({13.026207,  24.131974,  41.869426,   4.984373,  0.765769}), .c =   0.409294};
    constexpr FiveGaussian OH_alc      {.a = { 0.456221, 3.219608, 0.812773,  2.666928,  1.380927}, .b = s_to_q({21.503498,  13.397134,  34.547137,   5.826620,  0.412902}), .c =   0.463202};
    constexpr FiveGaussian NH          {.a = { 1.650531, 0.429639, 2.144736,  1.851894,  1.408921}, .b = s_to_q({10.603730,   6.987283,  29.939901,  10.573859,  0.611678}), .c =   0.510589};
    constexpr FiveGaussian NH2         {.a = { 1.904157, 1.942536, 2.435585,  0.730512,  1.379728}, .b = s_to_q({10.803702,  10.792421,  29.610479,   6.847755,  0.709687}), .c =   0.603738};
    constexpr FiveGaussian NH3_plus    {.a = { 1.882162, 1.933200, 2.465843,  0.927311,  1.190889}, .b = s_to_q({10.975157,  10.956008,  29.208572,   6.663555,  0.843650}), .c =   0.597322};
    constexpr FiveGaussian SH          {.a = { 0.570042, 6.337416, 1.641643,  5.398549,  1.527982}, .b = s_to_q({11.447986,   1.197657,  55.401032,  22.420955,  2.356552}), .c =   1.523944};
}

TEST_CASE("FormFactor::xray::raw: groups roughly reproduce the legacy Grudinin fits") {
    std::vector<std::pair<form_factor_t, FiveGaussian>> legacy = {
        {form_factor_t::CH,  grudinin::CH_sp3},
        {form_factor_t::CH2, grudinin::CH2_sp3},
        {form_factor_t::CH3, grudinin::CH3_sp3},
        {form_factor_t::OH,  grudinin::OH_alc},
        {form_factor_t::NH,  grudinin::NH},
        {form_factor_t::NH2, grudinin::NH2},
        {form_factor_t::NH3, grudinin::NH3_plus},
        {form_factor_t::SH,  grudinin::SH}
    };

    // the difference grows with q; CH3 is furthest off, at 7.4% at q = 1
    std::vector<std::pair<double, double>> tolerances = {{0.1, 0.005}, {0.3, 0.01}, {0.5, 0.025}, {0.7, 0.05}, {1.0, 0.1}};
    for (const auto& [type, coefficients] : legacy) {
        INFO("form factor " << to_string(type));
        xray::AtomicFormFactor old_ff(coefficients);
        const auto& new_ff = xray::raw::get(type);

        CHECK_THAT(new_ff.I0(), Catch::Matchers::WithinRel(old_ff.I0(), 1e-4));
        for (const auto& [q, tolerance] : tolerances) {
            INFO("q = " << q);
            CHECK_THAT(new_ff.evaluate(q), Catch::Matchers::WithinRel(old_ff.evaluate(q), tolerance));
        }
    }
}
