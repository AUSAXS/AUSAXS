#pragma once

#include <string_view>

namespace hist::detail::ff {
    // The form factor type of an atom. This is intended to be used as an index for best performance.
    enum class form_factor_t {
        HYDROGEN,
        CARBON,
        NITROGEN,
        OXYGEN,
        OTHER,
        COUNT
    };

    /**
     * @brief Get the number of unique form factors.
     */
    constexpr unsigned int get_count();

    /**
     * @brief Get the form factor type of an element from the Atom class. 
     */
    form_factor_t get_type(std::string_view element);

    constexpr double get_form_factor(form_factor_t ff_type1, form_factor_t ff_type2, double q);
}