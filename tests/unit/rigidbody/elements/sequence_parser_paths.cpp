// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>

#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/Rigidbody.h>
#include <data/Molecule.h>
#include <settings/All.h>
#include <io/ExistingFile.h>
#include <io/Folder.h>

#include <support/temp_file.h>

#include <filesystem>
#include <memory>
#include <string>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::sequencer;

// a scratch tree of its own, so the script can be given paths that only resolve relative to itself
struct PathResolutionFixture {
    PathResolutionFixture() {
        settings::general::verbose = false;
        settings::molecule::implicit_hydrogens = false;
        settings::grid::min_bins = 250;
        settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
        root = "temp/paths_" + test::detail::unique_tag();
        std::filesystem::create_directories(root + "/sub");
    }

    ~PathResolutionFixture() {
        settings::general::output = original_output;
        std::error_code ec;
        std::filesystem::remove_all(root, ec);
    }

    // place a copy of the test structure at <root>/sub/mol.pdb, reachable only through the config folder
    void plant_structure() const {
        std::filesystem::copy_file("tests/files/SASDJG5_single.pdb", root + "/sub/mol.pdb",
                                   std::filesystem::copy_options::overwrite_existing);
    }

    std::unique_ptr<Sequencer> parse(const std::string& contents) const {
        io::File config(io::Folder(root), "config", ".conf");
        config.create(contents);
        SequenceParser parser;
        return parser.parse_file(io::ExistingFile(config));
    }

    std::string root;
    std::string original_output = settings::general::output;
};

TEST_CASE_METHOD(PathResolutionFixture, "SequenceParser::LoadElement resolves paths relative to the script", "[files]") {
    SECTION("a structure in a subfolder next to the configuration file") {
        plant_structure();
        auto seq = parse(
            "load {\n"
            "    pdb sub/mol.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
        );
        REQUIRE(seq != nullptr);
        CHECK(seq->_get_rigidbody()->molecule.size_body() == 1);
    }

    SECTION("a structure sitting directly next to the configuration file") {
        std::filesystem::copy_file("tests/files/SASDJG5_single.pdb", root + "/mol.pdb",
                                   std::filesystem::copy_options::overwrite_existing);
        auto seq = parse(
            "load {\n"
            "    pdb mol.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
        );
        REQUIRE(seq != nullptr);
        CHECK(seq->_get_rigidbody()->molecule.size_body() == 1);
    }

    SECTION("a path that exists neither way is still an error") {
        CHECK_THROWS(parse(
            "load {\n"
            "    pdb sub/not_a_file.pdb\n"
            "    saxs tests/files/SASDJG5.dat\n"
            "}\n"
        ));
    }
}

TEST_CASE_METHOD(PathResolutionFixture, "SequenceParser::OutputFolderElement modes", "[files]") {
    plant_structure();
    static const std::string load =
        "load {\n"
        "    pdb sub/mol.pdb\n"
        "    saxs tests/files/SASDJG5.dat\n"
        "}\n";

    SECTION("absolute anchors the folder against the working directory") {
        REQUIRE_NOTHROW(parse(load + "output {\n    path results\n    mode absolute\n}\n"));
        CHECK(std::filesystem::path(settings::general::output).is_absolute());
    }

    SECTION("an already absolute path is left where it points") {
        auto target = std::filesystem::absolute(root).string() + "/results";
        REQUIRE_NOTHROW(parse(load + "output {\n    path " + target + "\n    mode absolute\n}\n"));
        CHECK(settings::general::output.starts_with(target));
    }

    SECTION("relative_config is prefixed with the script folder") {
        REQUIRE_NOTHROW(parse(load + "output {\n    path results\n    mode relative_config\n}\n"));
        CHECK(settings::general::output.starts_with(root));
    }

    SECTION("relative is left alone") {
        REQUIRE_NOTHROW(parse(load + "output results\n"));
        CHECK(settings::general::output == "results/");
    }

    SECTION("an unknown mode is still rejected") {
        CHECK_THROWS_AS(parse(load + "output {\n    path results\n    mode sideways\n}\n"), sequencer::except::parse_error);
    }
}
