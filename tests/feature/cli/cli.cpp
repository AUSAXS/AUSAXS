#include <catch2/catch_test_macros.hpp>

#include <constants/Version.h>
#include <io/File.h>
#include <support/cli_runner.h>
#include <support/temp_file.h>

#include <string>
#include <string_view>

using namespace ausaxs;

// These tests start the ausaxs executable as a subprocess, so they cover what libausaxs cannot see from the inside: argument handling, the exit code
// and what reaches the terminal.

TEST_CASE("cli: no tool prints the usage text and fails") {
    auto result = test::run_cli("");
    CHECK(result.exit_code == 1);
    CHECK(result.contains("Usage: ausaxs <tool>"));
    CHECK(result.contains("fit"));
    CHECK(result.contains("rigidbody"));
}

TEST_CASE("cli: an unknown tool is named back to the user") {
    auto result = test::run_cli("frobnicate");
    CHECK(result.exit_code == 1);
    CHECK(result.contains("Unknown tool: frobnicate"));
    CHECK(result.contains("Usage: ausaxs <tool>"));
}

TEST_CASE("cli: the top level answers --help and --version") {
    for (std::string_view flag : {"--help", "-h"}) {
        auto result = test::run_cli(flag);
        CHECK(result.exit_code == 0);
        CHECK(result.contains("Usage: ausaxs <tool>"));
        CHECK(!result.contains("Unknown tool"));
    }

    for (std::string_view flag : {"--version", "-v"}) {
        auto result = test::run_cli(flag);
        CHECK(result.exit_code == 0);
        CHECK(result.contains(constants::version));
        CHECK(!result.contains("Unknown tool"));
    }
}

TEST_CASE("cli: every tool answers --help successfully") {
    // a tool whose --help exits non-zero is unusable in a pipeline, and there is no other way to discover its options
    for (std::string_view tool : {"fit", "em", "rigidbody"}) {
        auto result = test::run_cli(std::string(tool) + " --help");
        CHECK(result.exit_code == 0);
        CHECK(!result.output.empty());
    }
}

TEST_CASE("cli: a missing input file is reported without a crash", "[files]") {
    auto result = test::run_cli("fit tests/files/no_such_structure.pdb tests/files/2epe.dat");
    CHECK(result.exit_code != 0);
    CHECK(result.exit_code < 128); // a signal would mean a crash
    CHECK(result.contains("no_such_structure.pdb"));
}

TEST_CASE("cli: a fit succeeds, exits zero and writes its results", "[files]") {
    test::TempFolder out;
    auto result = test::run_cli("fit tests/files/2epe.pdb tests/files/2epe.dat -o " + out.path() + "/");

    CHECK(result.exit_code == 0);
    CHECK(result.contains("chi2/dof"));

    // the tool appends the dataset stem to the output folder
    CHECK(io::File(out.path() + "/2epe/ausaxs.fit").exists());
    CHECK(io::File(out.path() + "/2epe/report.txt").exists());
}
