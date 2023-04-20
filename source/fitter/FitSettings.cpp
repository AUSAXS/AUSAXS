#include <fitter/FitSettings.h>

namespace settings::fit {
    settings::detail::SmartOption<bool> verbose(false, "fit-verbose");
    settings::detail::SmartOption<unsigned int> N(100, "N-points");
    settings::detail::SmartOption<unsigned int> max_iterations(100, "max-iterations");
}