/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <settings/CrystalSettings.h>
#include <settings/SettingsIORegistry.h>

unsigned int settings::crystal::h = 100;
unsigned int settings::crystal::k = 100;
unsigned int settings::crystal::l = 100;

double settings::crystal::max_q = 1e6; 
double settings::crystal::grid_expansion = 3;

double settings::crystal::reduced::basis_q = 3;

bool settings::crystal::detail::use_checkpointing = true;

namespace settings::crystal::io {
    settings::io::SettingSection grid_settings("Crystal", { 
        settings::io::create(h, "h"),
        settings::io::create(k, "k"),
        settings::io::create(l, "l"),
        settings::io::create(max_q, "max_q"),
        settings::io::create(grid_expansion, "grid_expansion"),
        settings::io::create(reduced::basis_q, "basis_q")
    });
}

namespace settings::crystal {
    MillerGenerationChoice miller_generation_strategy = MillerGenerationChoice::All;
}