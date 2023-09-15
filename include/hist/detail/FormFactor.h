#pragma once

#include <hist/detail/FormFactorType.h>
#include <hist/detail/FormFactor.h>
#include <utility/Constants.h>

#include <array>
#include <string_view>

namespace hist::detail {
    struct FormFactor {
        constexpr FormFactor(std::array<double, 5> a, std::array<double, 5> b, double c) : a(a), b(b), c(c) {}

        double evaluate(double q) const;

        /**
         * @brief Get the number of unique form factors.
         */
        static unsigned int get_count();

        /**
         * @brief Get the form factor type of an element from the Atom class. 
         */
        static form_factor_t get_type(std::string_view element);

        std::array<double, 5> a;
        std::array<double, 5> b;
        double c;
    };

    /**
     * This struct contains the form factors of the most common atomic elements encountered in SAXS. 
     */
    struct FormFactorStorage {
        static const FormFactor& get_form_factor(form_factor_t type);

        inline static const FormFactor hydrogen = FormFactor(constants::form_factor::hydrogen::a, constants::form_factor::hydrogen::b, constants::form_factor::hydrogen::c);
        inline static const FormFactor carbon   = FormFactor(  constants::form_factor::carbon::a,   constants::form_factor::carbon::b,   constants::form_factor::carbon::c);
        inline static const FormFactor nitrogen = FormFactor(constants::form_factor::nitrogen::a, constants::form_factor::nitrogen::b, constants::form_factor::nitrogen::c);
        inline static const FormFactor oxygen   = FormFactor(  constants::form_factor::oxygen::a,   constants::form_factor::oxygen::b,   constants::form_factor::oxygen::c);
        inline static const FormFactor other    = FormFactor(   constants::form_factor::other::a,    constants::form_factor::other::b,    constants::form_factor::other::c);
    };
}