#include <settings/GeneralSettings.h>
#include <settings/SettingsIORegistry.h>

#include <thread>

namespace settings::general {
    constexpr const char* const residue_folder = "temp/residues/";
    bool verbose = true;
    unsigned int threads = std::thread::hardware_concurrency();
    std::string output = "output/";
    bool keep_hydrogens = false;
    bool supplementary_plots = true;

    namespace detail {
        unsigned int job_size = 200; // The number of atoms to process in each job.
    };

    namespace io {
        settings::io::SettingSection general_settings("General", {
            settings::io::create(verbose, {"verbose", "v"}),
            settings::io::create(threads, {"threads", "t"}),
            settings::io::create(output, {"output", "o"}),
        });
    }
}