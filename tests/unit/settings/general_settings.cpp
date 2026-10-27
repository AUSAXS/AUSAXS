#include <catch2/catch_test_macros.hpp>

#include <settings/GeneralSettings.h>
#include <settings/SettingRef.h>
#include <settings/SettingsIO.h>

#include <string>

using namespace ausaxs;

// the residue cache must follow the cache setting, also when it is changed by name as a settings file does
TEST_CASE("settings::general::helper::residue_folder follows the cache setting") {
    std::string original = settings::general::cache;
    auto& stored = settings::io::detail::ISettingRef::get_stored_settings();

    stored.at("cache")->set({"/tmp/ausaxs_cache/"});
    CHECK(settings::general::helper::residue_folder() == "/tmp/ausaxs_cache/residues/");

    // a missing trailing separator is added
    stored.at("cache")->set({"/tmp/other_cache"});
    CHECK(settings::general::helper::residue_folder() == "/tmp/other_cache/residues/");

    settings::reset();
    CHECK(settings::general::helper::residue_folder() == original + "residues/");
}
