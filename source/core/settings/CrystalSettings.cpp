// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <settings/CrystalSettings.h>
#include <settings/SettingsIORegistry.h>

using namespace ausaxs;

unsigned int settings::crystal::h = 100;
unsigned int settings::crystal::k = 100;
unsigned int settings::crystal::l = 100;
double settings::crystal::max_q = 1e6; 
double settings::crystal::grid_expansion = 3;
double settings::crystal::reduced::basis_q = 3;
bool settings::crystal::detail::use_checkpointing = true;
ausaxs::settings::crystal::MillerGenerationChoice ausaxs::settings::crystal::miller_generation_strategy = ausaxs::settings::crystal::MillerGenerationChoice::All;

//? disabled for now since crystal calculations are not practically used anyway
// namespace ausaxs::settings::io {
//     settings::io::SettingSection crystal_section("Crystal", { 
//         settings::io::create(crystal::h, "h"),
//         settings::io::create(crystal::k, "k"),
//         settings::io::create(crystal::l, "l"),
//         settings::io::create(crystal::max_q, "max_q"),
//         settings::io::create(crystal::grid_expansion, "grid_expansion"),
//         settings::io::create(crystal::reduced::basis_q, "basis_q")
//     });
// }