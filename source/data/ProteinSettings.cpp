#include <data/ProteinSettings.h>

namespace settings::protein {
    settings::detail::SmartOption<bool> center(true, "center");
    settings::detail::SmartOption<bool> use_effective_charge(true, "use_effective_charge");
}