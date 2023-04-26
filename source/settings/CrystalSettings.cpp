#include <settings/CrystalSettings.h>
#include <settings/SettingsIORegistry.h>

namespace settings::crystal {
    unsigned int h = 100;
    unsigned int k = 100;
    unsigned int l = 100;

    double max_q = 1e6; 
    double grid_expansion = 3;

    double reduced::basis_q = 3;

    namespace io {
        settings::io::SettingSection grid_settings("Crystal", { 
            settings::io::create(h, "h"),
            settings::io::create(k, "k"),
            settings::io::create(l, "l"),
            settings::io::create(max_q, "max_q"),
            settings::io::create(grid_expansion, "grid_expansion"),
            settings::io::create(reduced::basis_q, "basis_q")
        });
    }
}

namespace settings::crystal {
    MillerGenerationChoice miller_generation_strategy = MillerGenerationChoice::All;
}