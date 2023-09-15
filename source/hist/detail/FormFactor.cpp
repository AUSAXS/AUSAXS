#include <hist/detail/FormFactor.h>

#include <cmath>

using namespace hist::detail;

const FormFactor& FormFactorStorage::get_form_factor(form_factor_t type) {
    switch (type) {
        case form_factor_t::HYDROGEN:
            return hydrogen;
        case form_factor_t::CARBON:
            return carbon;
        case form_factor_t::NITROGEN:
            return nitrogen;
        case form_factor_t::OXYGEN:
            return oxygen;
        case form_factor_t::OTHER:
            return other;
        default:
            throw std::runtime_error("Invalid form factor type");
    }
}

double FormFactor::evaluate(double q) const {
    double sum = 0;
    for (unsigned int i = 0; i < 5; ++i) {
        sum += a[i]*std::exp(-b[i]*q*q);
    }
    return sum + c;
}

form_factor_t FormFactor::get_type(std::string_view element) {
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

unsigned int FormFactor::get_count() {
    return static_cast<unsigned int>(form_factor_t::COUNT);
}