#include <form_factor/FormFactor.h>
#include <data/record/Atom.h>
#include <constants/ConstantsFwd.h>

#include <cmath>

using namespace form_factor;

constexpr const FormFactor& storage::get_form_factor(form_factor_t type) {
    switch (type) {
        case form_factor_t::H:
            return H;
        case form_factor_t::C:
            return C;
        case form_factor_t::N:
            return N;
        case form_factor_t::O:
            return O;
        case form_factor_t::S:
            return S;
        case form_factor_t::CH:
            return CH_sp3;
        case form_factor_t::CH2:
            return CH2_sp3;
        case form_factor_t::CH3:
            return CH3_sp3;
        case form_factor_t::NH:
            return NH;
        case form_factor_t::NH2:
            return NH2;
        case form_factor_t::OH:
            return OH_alc;
        case form_factor_t::SH:
            return SH;
        case form_factor_t::OTHER:
            return other;
        case form_factor_t::EXCLUDED_VOLUME:
            return excluded_volume;
        default:
            throw std::runtime_error("FormFactorStorage::get_form_factor: Invalid form factor type");
    }
}

double FormFactor::evaluate(double q) const {
    double sum = 0;
    for (unsigned int i = 0; i < 5; ++i) {
        sum += a[i]*std::exp(-b[i]*q*q);
    }
    return (sum + c)/f0;
}

constexpr form_factor_t form_factor::get_type(const data::record::Atom& atom) {
    switch(atom.element) {
        case constants::atom_t::H: return form_factor_t::H;
        case constants::atom_t::C: return form_factor_t::C;
        case constants::atom_t::N: return form_factor_t::N;
        case constants::atom_t::O: return form_factor_t::O;
        default: return form_factor_t::OTHER;
    }
}

constexpr unsigned int form_factor::get_count() {
    return static_cast<unsigned int>(form_factor_t::COUNT);
}

constexpr unsigned int form_factor::get_count_without_excluded_volume() {
    return static_cast<unsigned int>(form_factor_t::COUNT)-1;
}