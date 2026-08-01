// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>

#include <settings/SettingsIO.h>
#include <settings/GeneralSettings.h>
#include <io/File.h>

#include <support/temp_file.h>

using namespace ausaxs;

struct SettingsIOFixture {
    SettingsIOFixture() : original(settings::general::output) {}
    ~SettingsIOFixture() {settings::general::output = original;}

    // write the current settings to a scratch file, clear the output setting, and read them back
    static std::string round_trip(const std::string& value) {
        settings::general::output = value;
        test::TempFile file("ausaxs_settings_io_test", ".txt");
        settings::write(file);
        settings::general::output = "";
        settings::read(file);
        return settings::general::output;
    }

    std::string original;
};

// The output folder is a path, and a path may contain spaces. settings::read splits on whitespace,
// so settings::write has to quote such a value for it to survive the round trip.
TEST_CASE_METHOD(SettingsIOFixture, "settings: string values round-trip through write/read") {
    SECTION("no spaces") {
        CHECK(round_trip("output/saxs_fitter/") == "output/saxs_fitter/");
    }

    SECTION("spaces") {
        CHECK(round_trip("my output/saxs fitter/") == "my output/saxs fitter/");
    }

    SECTION("windows path with spaces") {
        CHECK(round_trip("C:\\Users\\John Doe\\out\\") == "C:\\Users\\John Doe\\out\\");
    }
}

TEST_CASE_METHOD(SettingsIOFixture, "settings: a hand-written quoted path is read as one value") {
    test::TempFile file("ausaxs_settings_io_test", ".txt", "output \"C:\\my folder\\out\\\"\n");
    settings::general::output = "";
    settings::read(file);
    CHECK(settings::general::output == "C:\\my folder\\out\\");
}
