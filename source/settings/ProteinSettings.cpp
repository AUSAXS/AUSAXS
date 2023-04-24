#include <settings/ProteinSettings.h>
#include <settings/SettingsIORegistry.h>

namespace settings::protein {
    bool center = true;
    bool use_effective_charge = true;

    namespace io {
        settings::io::SettingSection protein_settings("Protein", {
            settings::io::create(center, "center"),
            settings::io::create(use_effective_charge, "use_effective_charge")
        });
    }
}