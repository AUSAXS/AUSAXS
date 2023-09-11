#include <hist/detail/FormFactorType.h>

hist::detail::ff::form_factor_t hist::detail::ff::get_type(std::string_view element) {
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