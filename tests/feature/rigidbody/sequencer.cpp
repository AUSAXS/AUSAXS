#include <catch2/catch_test_macros.hpp>
#include <fitter/FitResult.h>  // IWYU pragma: keep
#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <settings/All.h>

#include <memory>
#include <string>

using namespace ausaxs;
using namespace ausaxs::rigidbody;

namespace {
    std::unique_ptr<sequencer::Sequencer> parse_sequence(const std::string& script) {
        return sequencer::SequenceParser().parse_text(script + (script.find("\nloop ") != std::string::npos ? "end\n" : ""));
    }

    void configure_settings() {
        settings::general::verbose = false;
        settings::grid::min_bins = 250;
        settings::molecule::implicit_hydrogens = false;
    }
}

TEST_CASE("Sequencer: parser basic run", "[files]") {
    configure_settings();
    auto seq = parse_sequence("load {\n    pdb tests/files/SASDJG5.pdb\n    saxs tests/files/SASDJG5.dat\n}\nloop 5\n    optimize_once\nend\n");
    auto result = seq->execute();
    REQUIRE(result != nullptr);
    CHECK(result->fval > 0);
}

TEST_CASE("Sequencer: parser split loading", "[files]") {
    configure_settings();
    REQUIRE_NOTHROW(parse_sequence("load {\n    pdb tests/files/LAR1-2.pdb\n    split 9 99\n    saxs tests/files/LAR1-2.pdb\n}\n"));
}

TEST_CASE("Sequencer: parser nested loops", "[files]") {
    configure_settings();
    auto seq = parse_sequence("load {\n    pdb tests/files/SASDJG5.pdb\n    saxs tests/files/SASDJG5.dat\n}\nloop 3\n    loop 2\n        optimize_once\n    end\nend\n");
    auto result = seq->execute();
    REQUIRE(result != nullptr);
    CHECK(result->fval > 0);
}

TEST_CASE("Sequencer: parser strategy configuration", "[files]") {
    configure_settings();
    auto seq = parse_sequence("load {\n    pdb tests/files/SASDJG5.pdb\n    saxs tests/files/SASDJG5.dat\n}\nloop 5\n    select random_body\n    transform rigid\n    optimize_once\nend\n");
    auto result = seq->execute();
    REQUIRE(result != nullptr);
    CHECK(result->fval > 0);
}

TEST_CASE("Sequencer: parser automatic constraints", "[files]") {
    configure_settings();
    auto seq = parse_sequence("load {\n    pdb tests/files/LAR1-2.pdb\n    split 9 99\n    saxs tests/files/LAR1-2.dat\n}\nautoconstrain backbone\nloop 5\n    optimize_once\nend\n");
    auto result = seq->execute();
    REQUIRE(result != nullptr);
    CHECK(result->fval > 0);
}
