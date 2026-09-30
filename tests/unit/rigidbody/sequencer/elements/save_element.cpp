#include <catch2/catch_test_macros.hpp>

#include <fitter/FitResult.h>  // IWYU pragma: keep
#include <io/Folder.h>
#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <settings/All.h>

#include <support/temp_file.h>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::sequencer;

struct SequencerElementsFixture {
    SequencerElementsFixture() {
        settings::general::verbose = false;
        settings::molecule::implicit_hydrogens = false;
        settings::grid::min_bins = 250;
    }
};

namespace {
    std::unique_ptr<Sequencer> parse_sequence(const std::string& script) {
        SequenceParser parser;
        return parser.parse_text(script + (script.find("\nloop ") != std::string::npos ? "end\n" : ""));
    }
}

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::SaveElement basic functionality") {
    SECTION("Save PDB file - verify no crash") {
        io::Folder out_dir("temp/ausaxs_test_output_" + test::detail::unique_tag());
        out_dir.create();
        std::string output_path = out_dir.path() + "/test_save.pdb";
        auto seq = parse_sequence(
            "load {\n"
            "    pdb tests/files/2epe.pdb\n"
            "    saxs tests/files/2epe.dat\n"
            "    split 40 80\n"
            "}\n"
            "loop 2\n"
            "    optimize_once\n"
            "    save " + output_path + "\n"
            "end\n"
        );

        REQUIRE(seq != nullptr);
        REQUIRE_NOTHROW(seq->execute());
    }
}
