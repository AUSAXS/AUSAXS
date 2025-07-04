// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <settings/MDSettings.h>
#include <settings//SettingsIORegistry.h>

using namespace ausaxs;

std::string settings::md::gmx_path = "gmx";
std::string settings::md::buffer_path = "";
std::string settings::md::water_model = "tip4p2005";
std::string settings::md::force_field = "amber99sb-ildn";
std::string settings::md::box_type = "dodecahedron";
std::string settings::md::cation = "Na";
std::string settings::md::anion = "Cl";
std::string settings::md::minimization_sim_location = "lucy";
std::string settings::md::thermalization_sim_location = "smaug";
std::string settings::md::production_sim_location = "smaug";

//? disabled; MD may be removed in the future
// namespace ausaxs::settings::io {
//     settings::io::SettingSection md_section("MD", {
//         settings::io::create(md::gmx_path, {"gmx_exe", "gmx_executable", "gmx"}),
//         settings::io::create(md::buffer_path, {"buffer_path", "buffer"}),
//         settings::io::create(md::water_model, {"water_model", "water"}),
//         settings::io::create(md::force_field, {"force_field", "forcefield", "ff"}),
//         settings::io::create(md::box_type, {"box_type", "boxtype", "box"}),
//         settings::io::create(md::cation, {"cation", "cation_type", "cationtype"}),
//         settings::io::create(md::anion, {"anion", "anion_type", "aniontype"}),
//         settings::io::create(md::minimization_sim_location, "minimization_sim_location"),
//         settings::io::create(md::thermalization_sim_location, "thermalization_sim_location"),
//         settings::io::create(md::production_sim_location, "production_sim_location")
//     });
// }