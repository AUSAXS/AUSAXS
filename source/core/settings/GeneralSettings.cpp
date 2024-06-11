/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <settings/GeneralSettings.h>
#include <settings/SettingsIORegistry.h>

#include <thread>

constexpr const char* const settings::general::residue_folder = "temp/residues/";
bool settings::general::verbose = true;
bool settings::general::warnings = true;
unsigned int settings::general::threads = std::thread::hardware_concurrency()-1;
std::string settings::general::output = "output/";
bool settings::general::keep_hydrogens = false;
bool settings::general::supplementary_plots = true;

namespace settings::general::detail {
    unsigned int job_size = 800; // The number of atoms to process in each job.
};

namespace settings::general::io {
    settings::io::SettingSection general_settings("General", {
        settings::io::create(verbose, {"verbose", "v"}),
        settings::io::create(warnings, {"warnings", "w"}),
        settings::io::create(threads, {"threads", "t"}),
        settings::io::create(output, {"output", "o"}),
    });
}