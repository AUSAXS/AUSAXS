// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <settings/GeneralSettings.h>
#include <settings/SettingsIORegistry.h>
#include <utility/StringUtils.h>

#include <thread>
#include <algorithm>

using namespace ausaxs;

bool settings::general::verbose = true;
bool settings::general::warnings = true;
bool settings::general::offline = false;
unsigned int settings::general::threads = std::thread::hardware_concurrency()-1;
std::string settings::general::output = "output/";
bool settings::general::keep_hydrogens = false;
bool settings::general::supplementary_plots = true;
bool settings::general::generate_plots = true;
settings::general::QUnit settings::general::input_q_unit = settings::general::QUnit::A;

std::string settings::general::cache = [] () {
    const char* env_p = nullptr;

    #ifdef _WIN32
        env_p = std::getenv("LOCALAPPDATA");
        if (env_p) {
            return std::string(env_p) + "/ausaxs/";
        }
    #elif defined(__APPLE__)
        env_p = std::getenv("HOME");
        if (env_p) {
            return std::string(env_p) + "/Library/Caches/ausaxs/";
        }
    #else
        env_p = std::getenv("XDG_CACHE_HOME");
        if (env_p) {
            return std::string(env_p) + "/ausaxs/";
        }

        env_p = std::getenv("HOME");
        if (env_p) {
            return std::string(env_p) + "/.cache/ausaxs/";
        }
    #endif
    return output + "/temp/";
}();
std::string settings::general::residue_folder = cache + "residues/";

unsigned int ausaxs::settings::general::detail::get_job_size(unsigned int n) {
    constexpr unsigned int jobs_per_thread = 8; // aim for this many jobs per thread
    constexpr unsigned int min_job_size = 64; // but never go below this to avoid excessive overhead from too many small jobs
    constexpr unsigned int max_job_size = 512;

    unsigned int threads = std::max(1u, settings::general::threads);
    return std::clamp(n/(jobs_per_thread*threads), min_job_size, max_job_size);
}

namespace ausaxs::settings::io {
    settings::io::SettingSection general_section("General", {
        settings::io::create(general::verbose, {"verbose", "v"}),
        settings::io::create(general::warnings, {"warnings", "w"}),
        settings::io::create(general::threads, {"threads", "t"}),
        settings::io::create(general::output, {"output", "o"}),
        settings::io::create(general::cache, {"cache"}),
        settings::io::create(general::keep_hydrogens, {"keep_hydrogens"}),
        settings::io::create(general::supplementary_plots, {"supplementary_plots"}),
        settings::io::create(general::input_q_unit, {"unit"}),
        settings::io::create(general::offline, {"offline"})
    });
}

template<> std::string settings::io::detail::SettingRef<settings::general::QUnit>::get() const {
    switch (settingref) {
        case settings::general::QUnit::A: return "default_to_A";
        case settings::general::QUnit::NM: return "default_to_nm";
        case settings::general::QUnit::USER_A: return "A";
        case settings::general::QUnit::USER_NM: return "nm";
        default: return std::to_string(static_cast<int>(settingref));
    }
}

template<> void settings::io::detail::SettingRef<settings::general::QUnit>::set(const std::vector<std::string>& val) {
    auto str = utility::to_lowercase(val[0]);
    if (str == "a" || str == "å") {settingref = settings::general::QUnit::USER_A;}
    else if (str == "nm") {settingref = settings::general::QUnit::USER_NM;}
    else if (str == "default_to_a") {settingref = settings::general::QUnit::A;}
    else if (str == "default_to_nm") {settingref = settings::general::QUnit::NM;}
    else {throw ausaxs::except::invalid_argument("settings: Unknown value for q-unit: \"" + val[0] + "\"");}
}

template<> std::string settings::io::detail::SettingRef<settings::general::QUnit>::type() const {return "string";}

bool ausaxs::settings::general::helper::is_angstroms(QUnit u) {
    // QUnit::USER_A follows from QUnit::A due to the bit manipulation
    return static_cast<char>(u) & static_cast<char>(QUnit::A);
}

bool ausaxs::settings::general::helper::is_nanometers(QUnit u) {
    // QUnit::USER_NM follows from QUnit::NM due to the bit manipulation
    return static_cast<char>(u) & static_cast<char>(QUnit::NM);
}

bool ausaxs::settings::general::helper::is_user_defined(QUnit u) {
    return static_cast<char>(u) & (1 << 3);
}