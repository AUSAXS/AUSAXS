#include <settings/MoleculeSettings.h>
#include <settings/SettingsIORegistry.h>

bool settings::molecule::center = true;
bool settings::molecule::use_effective_charge = true;

namespace settings::molecule::io {
    settings::io::SettingSection molecule_settings("Molecule", {
        settings::io::create(center, "center"),
        settings::io::create(use_effective_charge, "use_effective_charge")
    });
}