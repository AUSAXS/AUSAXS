#include <settings/MoleculeSettings.h>
#include <settings/SettingsIORegistry.h>

namespace settings::protein {
    bool center = true;
    bool use_effective_charge = true;

    namespace io {
        settings::io::SettingSection molecule_settings("Molecule", {
            settings::io::create(center, "center"),
            settings::io::create(use_effective_charge, "use_effective_charge")
        });
    }
}