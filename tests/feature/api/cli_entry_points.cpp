// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>

#include <api/cli/cli_em_fitter.h>
#include <api/cli/cli_rigidbody.h>
#include <api/cli/cli_saxs_fitter.h>
#include <settings/All.h>

#include <vector>

using namespace ausaxs;

// pyAUSAXS and its GUI call the exported CLI functions in-process, so every outcome must be a return code:
// an exit() or an escaping exception here would take the test binary down with it.
namespace {
    template<typename F>
    int call(F&& f, std::vector<const char*> argv) {
        return f(static_cast<int>(argv.size()), argv.data());
    }
}

TEST_CASE("cli: exported entry points return instead of exiting or throwing") {
    settings::general::verbose = false;

    SECTION("--version returns 0") {
        CHECK(call(cli_saxs_fitter, {"fit", "--version"}) == 0);
        CHECK(call(cli_em_fitter, {"em", "--version"}) == 0);
        CHECK(call(cli_rigidbody, {"rigidbody", "--version"}) == 0);
    }

    SECTION("a missing required input returns 1") {
        CHECK(call(cli_saxs_fitter, {"fit"}) == 1);
        CHECK(call(cli_em_fitter, {"em"}) == 1);
        CHECK(call(cli_rigidbody, {"rigidbody"}) == 1);
    }

    SECTION("an option callback that throws while parsing returns 1") {
        CHECK(call(cli_saxs_fitter, {"fit", "tests/files/2epe.pdb", "tests/files/2epe.dat", "exv", "--model", "Foo"}) == 1);
    }
}
