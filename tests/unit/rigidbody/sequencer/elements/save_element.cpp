#include <catch2/catch_test_macros.hpp>

#include <io/File.h>
#include <io/Folder.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/detail/parse_error.h>
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

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::SaveElement checks the format when parsed") {
    std::string load =
        "load {\n"
        "    pdb tests/files/2epe.pdb\n"
        "    saxs tests/files/2epe.dat\n"
        "    split 40 80\n"
        "}\n";
    for (std::string ok : {"a.pdb", "a_%.fit", "a.txt", "a.png", "a.xyz"}) {
        INFO(ok);
        CHECK_NOTHROW(parse_sequence(load + "save " + ok + "\n"));
    }

    // these used to parse, run the whole optimization, and only then throw when the save was reached
    for (std::string bad : {"final.dat", "final", "final.pdb.gz"}) {
        INFO(bad);
        CHECK_THROWS_AS(parse_sequence(load + "save " + bad + "\n"), rigidbody::sequencer::except::parse_error);
    }
}

TEST_CASE_METHOD(SequencerElementsFixture, "SequencerElements::SaveElement numbers each path template on its own") {
    io::Folder out_dir("temp/ausaxs_test_output_" + test::detail::unique_tag());
    out_dir.create();
    settings::general::output = out_dir.path() + "/"; // save paths are resolved relative to the output folder
    auto seq = parse_sequence(
        "load {\n"
        "    pdb tests/files/2epe.pdb\n"
        "    saxs tests/files/2epe.dat\n"
        "    split 40 80\n"
        "}\n"
        "save first_%.pdb\n"
        "loop 2\n"
        "    optimize_once\n"
        "    save every_%.pdb\n"
        "    save first_%.pdb\n"
        "end\n"
    );
    REQUIRE_NOTHROW(seq->execute());

    // the counters used to be shared per file type, so the second series started where the first left off
    auto exists = [&] (const std::string& name) {return io::File(out_dir.path() + "/" + name).exists();};
    CHECK(exists("every_1.pdb"));
    CHECK(exists("every_2.pdb"));
    CHECK_FALSE(exists("every_3.pdb"));

    // two elements saving to the same template continue one series
    CHECK(exists("first_1.pdb"));
    CHECK(exists("first_2.pdb"));
    CHECK(exists("first_3.pdb"));
}
