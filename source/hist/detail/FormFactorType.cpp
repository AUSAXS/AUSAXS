#include <hist/detail/FormFactorType.h>
#include <utility/Constants.h>

#include <cmath>

using namespace hist::detail::ff;

form_factor_t hist::detail::ff::get_type(std::string_view element) {
    if (element == "H") {
        return form_factor_t::HYDROGEN;
    } else if (element == "C") {
        return form_factor_t::CARBON;
    } else if (element == "N") {
        return form_factor_t::NITROGEN;
    } else if (element == "O") {
        return form_factor_t::OXYGEN;
    } else {
        return form_factor_t::OTHER;
    }
}

constexpr unsigned int hist::detail::ff::get_count() {
    return static_cast<unsigned int>(form_factor_t::COUNT);
}

EvalFormFactor::EvalFormFactor(form_factor_t ff1, form_factor_t ff2) {
    auto get_sigma = [](form_factor_t ff_type) {
        switch (ff_type) {
            case form_factor_t::HYDROGEN:
                return constants::form_factor::hydrogen;
            case form_factor_t::CARBON:
                return constants::form_factor::carbon;
            case form_factor_t::NITROGEN:
                return constants::form_factor::nitrogen;
            case form_factor_t::OXYGEN:
                return constants::form_factor::oxygen;
            case form_factor_t::OTHER:
                return constants::form_factor::other;
            default:
                throw std::runtime_error("EvalFormFactor::EvalFormFactor: Invalid form factor type.");
        }
    };
    factor = get_sigma(ff1)*get_sigma(ff2)/3;
}

double EvalFormFactor::operator()(double q) const {
    return std::exp(-factor*q*q);
}

constexpr EvalFormFactor hist::detail::ff::get_form_factor(int ff_type1, int ff_type2) {
    return EvalFormFactor(static_cast<form_factor_t>(ff_type1), static_cast<form_factor_t>(ff_type2));
}