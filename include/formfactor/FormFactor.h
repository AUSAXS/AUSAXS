#pragma once

#include <formfactor/FormFactorType.h>
#include <formfactor/FormFactor.h>
#include <formfactor/FormFactorTable.h>

#include <array>
#include <string_view>

class Atom;
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

            /**
             * @brief Get the number of unique form factors.
             * 
             * This can be used to iterate over all form factors.
             */
            static unsigned int get_count();

            /**
             * @brief Get the number of unique form factors, excluding the excluded volume form factor.
             * 
             * This can be used to iterate over all form factors except the excluded volume form factor.
             */
            static unsigned int get_count_without_excluded_volume();

            /**
             * @brief Get the form factor type of an Atom. 
             */
            static form_factor_t get_type(const Atom& atom);

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
     * This struct contains the form factors of the most common atomic elements encountered in SAXS. 
     */
    namespace storage {
        constexpr const FormFactor& get_form_factor(form_factor_t type);
        inline constexpr FormFactor hydrogen        = FormFactor(        constants::form_factor::hydrogen::a,         constants::form_factor::hydrogen::b,         constants::form_factor::hydrogen::c);
        inline constexpr FormFactor carbon          = FormFactor(  constants::form_factor::neutral_carbon::a,   constants::form_factor::neutral_carbon::b,   constants::form_factor::neutral_carbon::c);
        inline constexpr FormFactor nitrogen        = FormFactor(constants::form_factor::neutral_nitrogen::a, constants::form_factor::neutral_nitrogen::b, constants::form_factor::neutral_nitrogen::c);
        inline constexpr FormFactor oxygen          = FormFactor(  constants::form_factor::neutral_oxygen::a,   constants::form_factor::neutral_oxygen::b,   constants::form_factor::neutral_oxygen::c);
        inline constexpr FormFactor other           = FormFactor(           constants::form_factor::other::a,            constants::form_factor::other::b,            constants::form_factor::other::c);
        inline constexpr FormFactor excluded_volume = FormFactor( constants::form_factor::excluded_volume::a,  constants::form_factor::excluded_volume::b,  constants::form_factor::excluded_volume::c);
    };
}