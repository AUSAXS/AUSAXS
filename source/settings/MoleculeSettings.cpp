#include <settings/MoleculeSettings.h>
#include <settings/SettingsIORegistry.h>

bool settings::molecule::center = true;
bool settings::molecule::implicit_hydrogens = true;
bool settings::molecule::use_effective_charge = true;

#if DEBUG
    bool settings::molecule::throw_on_unknown_atom = true;
#else
    bool settings::molecule::throw_on_unknown_atom = false;
#endif

namespace settings::molecule::io {
    settings::io::SettingSection molecule_settings("Molecule", {
        settings::io::create(center, "center"),
        settings::io::create(use_effective_charge, "use_effective_charge"),
        settings::io::create(throw_on_unknown_atom, "throw_on_unknown_atom"),
        settings::io::create(implicit_hydrogens, "implicit_hydrogens")
    });
}