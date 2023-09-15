#pragma once

namespace hist::detail {
    // The form factor type of an atom. This is intended to be used as an index for best performance.
    enum class form_factor_t {
        HYDROGEN,
        CARBON,
        NITROGEN,
        OXYGEN,
        OTHER,
        COUNT
    };
}