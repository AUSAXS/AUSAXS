#include <catch2/catch_test_macros.hpp>

#include <fitter/FitResult.h>  // IWYU pragma: keep
#include <io/ExistingFile.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <settings/All.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody::sequencer;

TEST_CASE("SequenceParser: parse minimal config", "[files]") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::grid::min_bins = 250;

    std::string config =
        "load {\n"
        "    pdb tests/files/SASDJG5.pdb\n"
        "    saxs tests/files/SASDJG5.dat\n"
        "    split chain\n"
        "}\n"
        "loop 3\n"
        "    optimize_once\n"
        "    end\n"
        "end\n";

    SequenceParser parser;
    auto sequencer = parser.parse_text(config);
    REQUIRE(sequencer != nullptr);

    auto result = sequencer->execute();
    REQUIRE(result != nullptr);
    CHECK(result->fval > 0);
}

TEST_CASE("SequenceParser: parse normal config with output folder", "[files]") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::grid::min_bins = 250;

    SequenceParser parser;
    auto sequencer = parser.parse_file("tests/files/rigidbody/normal.conf");
    REQUIRE(sequencer != nullptr);

    auto result = sequencer->execute();
    REQUIRE(result != nullptr);
    CHECK(result->fval > 0);
}

TEST_CASE("SequenceParser: parse symmetry config", "[files]") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::grid::min_bins = 250;

    SequenceParser parser;
    auto sequencer = parser.parse_file("tests/files/rigidbody/symmetry.conf");
    REQUIRE(sequencer != nullptr);

    auto result = sequencer->execute();
    REQUIRE(result != nullptr);
    CHECK(result->fval > 0);
}

TEST_CASE("SequenceParser: a bare loop takes the nearest preceding parameter count", "[files]") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::grid::min_bins = 250;

    std::string config =
        "load {\n"
        "    pdb tests/files/LAR1-2.pdb\n"
        "    saxs tests/files/LAR1-2.dat\n"
        "    split 9, 99\n"
        "}\n"
        "parameter_strategy {\n"
        "    iterations 2\n"
        "    translate 1\n"
        "    rotate 1\n"
        "}\n"
        "loop\n"
        "    optimize_once\n"
        "    end\n"
        "end\n"
        "parameter_strategy {\n"
        "    iterations 5\n"
        "    translate 1\n"
        "    rotate 1\n"
        "}\n"
        "loop\n"
        "    optimize_once\n"
        "    end\n"
        "end\n";

    auto sequencer = SequenceParser().parse_text(config);
    std::vector<int> counts;
    for (auto& e : sequencer->_get_elements()) {
        if (auto* loop = dynamic_cast<LoopElement*>(e.get())) {counts.push_back(loop->_get_loop_iterations());}
    }
    CHECK(counts == std::vector<int>{2, 5});
}
