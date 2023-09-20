#pragma once

#include <hist/detail/FormFactorType.h>
#include <hist/detail/FormFactor.h>
#include <utility/Constants.h>

#include <array>
#include <string_view>

class Atom;
namespace hist::detail {
    class FormFactor {
        public:
            constexpr FormFactor(std::array<double, 5> a, std::array<double, 5> b, double c) : a(a), b(b), c(c) {
                normalization_factor = evaluate(0);
            }

            double evaluate(double q) const;

            /**
             * @brief Get the number of unique form factors.
             */
            static unsigned int get_count();

            /**
             * @brief Get the form factor type of an Atom. 
             */
            static form_factor_t get_type(const Atom& atom);

            std::array<double, 5> a;
            std::array<double, 5> b;
            double c;

        private: 
            double normalization_factor = 1;
    };

    /**
     * This struct contains the form factors of the most common atomic elements encountered in SAXS. 
     */
    struct FormFactorStorage {
        static const FormFactor& get_form_factor(form_factor_t type);

        inline static const FormFactor hydrogen        = FormFactor(        constants::form_factor::hydrogen::a,         constants::form_factor::hydrogen::b,         constants::form_factor::hydrogen::c);
        inline static const FormFactor carbon          = FormFactor(  constants::form_factor::neutral_carbon::a,   constants::form_factor::neutral_carbon::b,   constants::form_factor::neutral_carbon::c);
        inline static const FormFactor nitrogen        = FormFactor(constants::form_factor::neutral_nitrogen::a, constants::form_factor::neutral_nitrogen::b, constants::form_factor::neutral_nitrogen::c);
        inline static const FormFactor oxygen          = FormFactor(  constants::form_factor::neutral_oxygen::a,   constants::form_factor::neutral_oxygen::b,   constants::form_factor::neutral_oxygen::c);
        inline static const FormFactor other           = FormFactor(           constants::form_factor::other::a,            constants::form_factor::other::b,            constants::form_factor::other::c);
        inline static const FormFactor excluded_volume = FormFactor( constants::form_factor::excluded_volume::a,  constants::form_factor::excluded_volume::b,  constants::form_factor::excluded_volume::c);
    };
}