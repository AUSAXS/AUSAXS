#include <hist/detail/FormFactor.h>
#include <data/Atom.h>

#include <cmath>

using namespace hist::detail;

const FormFactor& FormFactorStorage::get_form_factor(form_factor_t type) {
    switch (type) {
        case form_factor_t::NEUTRAL_HYDROGEN:
            return hydrogen;
        case form_factor_t::NEUTRAL_CARBON:
            return carbon;
        case form_factor_t::NEUTRAL_NITROGEN:
            return nitrogen;
        case form_factor_t::NEUTRAL_OXYGEN:
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
    return sum + c - 5.95*std::exp(-1.62*q*q/2);
}

form_factor_t FormFactor::get_type(const Atom& atom) {
    auto element = atom.element;
    if (element == "H") {
        return form_factor_t::NEUTRAL_HYDROGEN;
    } else if (element == "C") {
        return form_factor_t::NEUTRAL_CARBON;
    } else if (element == "N") {
        return form_factor_t::NEUTRAL_NITROGEN;
    } else if (atom.is_water()) {
        return form_factor_t::NEUTRAL_OXYGEN;
    } else {
        return form_factor_t::OTHER;
    }
}

unsigned int FormFactor::get_count() {
    return static_cast<unsigned int>(form_factor_t::COUNT);
}