#pragma once

#include <form_factor/FormFactorType.h>
#include <form_factor/FormFactor.h>
#include <form_factor/FormFactorTable.h>
#include <data/DataFwd.h>

#include <array>
#include <string_view>

namespace form_factor {
    class FormFactor {
        public:
            constexpr FormFactor(std::array<double, 5> a, std::array<double, 5> b, double c) : a(a), b(b), c(c) {
                initialize();
            }

            /**
             * @brief Evaluate the form factor at a given q value.
             *        The form factor is normalized to 1 at q = 0.
             */
            double evaluate(double q) const;

            std::array<double, 5> a;
            std::array<double, 5> b;
            double c;

        private: 
            double f0 = 1;

            constexpr void initialize() {
                f0 = a[0] + a[1] + a[2] + a[3] + a[4] + c;
            }
    };


    /**
     * @brief Get the number of unique form factors.
     * 
     * This can be used to iterate over all form factors.
     */
    constexpr unsigned int get_count();

    /**
     * @brief Get the number of unique form factors, excluding the excluded volume form factor.
     * 
     * This can be used to iterate over all form factors except the excluded volume form factor.
     */
    constexpr unsigned int get_count_without_excluded_volume();

    /**
     * @brief Get the form factor type of an Atom. 
     */
    constexpr form_factor_t get_type(const data::record::Atom& atom);

    /**
     * This struct contains the form factors of the most common atomic elements encountered in SAXS. 
     */
    namespace storage {
        constexpr const FormFactor& get_form_factor(form_factor_t type);

        // atomic
        inline constexpr FormFactor H               = FormFactor(               constants::form_factor::H::a,               constants::form_factor::H::b,               constants::form_factor::H::c);
        inline constexpr FormFactor C               = FormFactor(               constants::form_factor::C::a,               constants::form_factor::C::b,               constants::form_factor::C::c);
        inline constexpr FormFactor N               = FormFactor(               constants::form_factor::N::a,               constants::form_factor::N::b,               constants::form_factor::N::c);
        inline constexpr FormFactor O               = FormFactor(               constants::form_factor::O::a,               constants::form_factor::O::b,               constants::form_factor::O::c);
        inline constexpr FormFactor S               = FormFactor(               constants::form_factor::S::a,               constants::form_factor::S::b,               constants::form_factor::S::c);

        // atomic groups
        inline constexpr FormFactor CH_sp3          = FormFactor(          constants::form_factor::CH_sp3::a,          constants::form_factor::CH_sp3::b,          constants::form_factor::CH_sp3::c);
        inline constexpr FormFactor CH2_sp3         = FormFactor(         constants::form_factor::CH2_sp3::a,         constants::form_factor::CH2_sp3::b,         constants::form_factor::CH2_sp3::c);
        inline constexpr FormFactor CH3_sp3         = FormFactor(         constants::form_factor::CH3_sp3::a,         constants::form_factor::CH3_sp3::b,         constants::form_factor::CH3_sp3::c);
        inline constexpr FormFactor CH_sp2          = FormFactor(          constants::form_factor::CH_sp2::a,          constants::form_factor::CH_sp2::b,          constants::form_factor::CH_sp2::c);
        inline constexpr FormFactor CH_arom         = FormFactor(         constants::form_factor::CH_arom::a,         constants::form_factor::CH_arom::b,         constants::form_factor::CH_arom::c);
        inline constexpr FormFactor OH_alc          = FormFactor(          constants::form_factor::OH_alc::a,          constants::form_factor::OH_alc::b,          constants::form_factor::OH_alc::c);
        inline constexpr FormFactor OH_acid         = FormFactor(         constants::form_factor::OH_acid::a,         constants::form_factor::OH_acid::b,         constants::form_factor::OH_acid::c);
        inline constexpr FormFactor O_res           = FormFactor(           constants::form_factor::O_res::a,           constants::form_factor::O_res::b,           constants::form_factor::O_res::c);
        inline constexpr FormFactor NH              = FormFactor(              constants::form_factor::NH::a,              constants::form_factor::NH::b,              constants::form_factor::NH::c);
        inline constexpr FormFactor NH2             = FormFactor(             constants::form_factor::NH2::a,             constants::form_factor::NH2::b,             constants::form_factor::NH2::c);
        inline constexpr FormFactor NH_plus         = FormFactor(         constants::form_factor::NH_plus::a,         constants::form_factor::NH_plus::b,         constants::form_factor::NH_plus::c);
        inline constexpr FormFactor NH2_plus        = FormFactor(        constants::form_factor::NH2_plus::a,        constants::form_factor::NH2_plus::b,        constants::form_factor::NH2_plus::c);
        inline constexpr FormFactor NH3_plus        = FormFactor(        constants::form_factor::NH3_plus::a,        constants::form_factor::NH3_plus::b,        constants::form_factor::NH3_plus::c);
        inline constexpr FormFactor NH_guanine      = FormFactor(      constants::form_factor::NH_guanine::a,      constants::form_factor::NH_guanine::b,      constants::form_factor::NH_guanine::c);
        inline constexpr FormFactor NH2_guanine     = FormFactor(     constants::form_factor::NH2_guanine::a,     constants::form_factor::NH2_guanine::b,     constants::form_factor::NH2_guanine::c);
        inline constexpr FormFactor SH              = FormFactor(              constants::form_factor::SH::a,              constants::form_factor::SH::b,              constants::form_factor::SH::c);

        // average excluded volume
        inline constexpr FormFactor excluded_volume = FormFactor( constants::form_factor::excluded_volume::a, constants::form_factor::excluded_volume::b, constants::form_factor::excluded_volume::c);

        // all others; this is just the form factor of argon
        inline constexpr FormFactor other           = FormFactor(           constants::form_factor::other::a,           constants::form_factor::other::b,           constants::form_factor::other::c);
    };
}