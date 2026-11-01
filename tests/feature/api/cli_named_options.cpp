// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>

#include <api/cli/cli_saxs_fitter.h>
#include <settings/All.h>

#include <string>
#include <vector>

using namespace ausaxs;

// a settings file may store an enum by its index, but on the command line a number used to be cast into the enum unchecked
TEST_CASE("cli: enum options refuse a number") {
    settings::general::verbose = false;
    auto exv_method = settings::exv::exv_method.value;
    auto exv_set = settings::exv::exv_set.value;
    auto hydration_strategy = settings::hydrate::hydration_strategy.value;

    for (const std::vector<std::string>& option : std::vector<std::vector<std::string>>{
        {"exv", "--model", "2"}, {"exv", "--model", "99"}, {"exv", "--table", "1"}, {"solv", "--model", "0"}
    }) {
        std::vector<const char*> argv = {"fit", "tests/files/2epe.pdb", "tests/files/2epe.dat"};
        for (const auto& arg : option) {argv.push_back(arg.c_str());}
        INFO(option[0] + " " + option[1] + " " + option[2]);
        CHECK(cli_saxs_fitter(static_cast<int>(argv.size()), argv.data()) != 0);
    }

    CHECK(settings::exv::exv_method == exv_method);
    CHECK(settings::exv::exv_set == exv_set);
    CHECK(settings::hydrate::hydration_strategy == hydration_strategy);
}
