#include <catch2/catch_test_macros.hpp>

#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <fitter/FitResult.h>
#include <settings/All.h>
#include <io/ExistingFile.h>

#include <fstream>

using namespace ausaxs;
using namespace ausaxs::rigidbody::sequencer;

TEST_CASE("SequenceParser: parse minimal config", "[files]") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::grid::min_bins = 500;

    // Write a minimal config to a temporary location
    std::string config_path = "/tmp/ausaxs_test_minimal.conf";
    {
        std::ofstream f(config_path);
        f << "load {\n"
          << "    pdb tests/files/SASDJG5.pdb\n"
          << "    saxs tests/files/SASDJG5.dat\n"
          << "    split chain\n"
          << "}\n"
          << "loop 3\n"
          << "    optimize_once\n"
          << "    end\n"
          << "end\n";
    }

    SequenceParser parser;
    auto sequencer = parser.parse(io::ExistingFile(config_path));
    REQUIRE(sequencer != nullptr);

    auto result = sequencer->execute();
    REQUIRE(result != nullptr);
    CHECK(result->fval > 0);
}
