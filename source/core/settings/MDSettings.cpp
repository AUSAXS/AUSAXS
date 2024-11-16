/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

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

namespace ausaxs::settings::md::io {
    settings::io::SettingSection md_section("MD", {
        settings::io::create(gmx_path, {"gmx_exe", "gmx_executable", "gmx"}),
        settings::io::create(buffer_path, {"buffer_path", "buffer"}),
        settings::io::create(water_model, {"water_model", "water"}),
        settings::io::create(force_field, {"force_field", "forcefield", "ff"}),
        settings::io::create(box_type, {"box_type", "boxtype", "box"}),
        settings::io::create(cation, {"cation", "cation_type", "cationtype"}),
        settings::io::create(anion, {"anion", "anion_type", "aniontype"}),
        settings::io::create(minimization_sim_location, "minimization_sim_location"),
        settings::io::create(thermalization_sim_location, "thermalization_sim_location"),
        settings::io::create(production_sim_location, "production_sim_location")
    });
}