#pragma once

#include <utility/SmartOption.h>

namespace settings {
    namespace protein {
        extern settings::detail::SmartOption<bool> center;               // Decides if the structure will be centered at origo.
        extern settings::detail::SmartOption<bool> use_effective_charge; // Decides whether the charge of the displaced water will be included.
    }
}