// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>

#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/elements/setup/RelativeHydrationElement.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/Rigidbody.h>
#include <data/Molecule.h>
#include <settings/All.h>
#include <io/ExistingFile.h>

#include <fstream>
#include <vector>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::sequencer;

struct SequenceParserRelativeHydrationFixture {
    SequenceParserRelativeHydrationFixture() {
        settings::general::verbose = false;
        settings::molecule::implicit_hydrogens = false;
        settings::grid::min_bins = 100;
        settings::hydrate::hydration_strategy = settings::hydrate::HydrationStrategy::NoStrategy;
        settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    }

    std::unique_ptr<Sequencer> parse(const std::string& content) {
        static int counter = 0;
        std::string path = "/tmp/ausaxs_seq_rel_hydration_test_" + std::to_string(counter++) + ".conf";
        std::ofstream f(path);
        f << content;
        f.close();
        SequenceParser parser;
        return parser.parse_file(path);
    }

    // 2epe split into three bodies at load time, plus whatever script lines follow
    std::unique_ptr<Sequencer> build(const std::string& extra = "") {
        return parse(
            "load {\n"
            "    pdb tests/files/2epe.pdb\n"
            "    saxs tests/files/2epe.dat\n"
            "    split 40 80\n"
            "}\n"
            + extra
        );
    }
};

TEST_CASE_METHOD(SequenceParserRelativeHydrationFixture, "RelativeHydrationElement: ratios are laid out over the whole body set", "[files]") {
    // BodyCounterCulling::cull indexes body_ratios by body index for every body in the molecule, so the vector must be exactly that long. It used to be built
    // by appending encoded name-registry indices onto the caller's ratio list, leaving it double length with indices standing in for weights.
    constexpr double normal = 1.0, high = 1.5, maximum = 1.75;

    SECTION("an unnamed body keeps the normal weight") {
        auto seq = build();
        REQUIRE(seq->_get_rigidbody()->molecule.size_body() == 3);

        RelativeHydrationElement element(seq.get(), {"b2"}, {maximum});
        CHECK(element._get_ratios() == std::vector<double>{normal, maximum, normal});
    }

    SECTION("several named bodies each land at their own index") {
        auto seq = build();
        RelativeHydrationElement element(seq.get(), {"b3", "b1"}, {high, maximum});
        CHECK(element._get_ratios() == std::vector<double>{maximum, normal, high});
    }

    SECTION("naming every body leaves no default weights") {
        auto seq = build();
        RelativeHydrationElement element(seq.get(), {"b1", "b2", "b3"}, {maximum, high, normal});
        CHECK(element._get_ratios() == std::vector<double>{maximum, high, normal});
    }

    SECTION("a renamed body is addressable by its alias") {
        auto seq = build("rename b2 core\n");
        RelativeHydrationElement element(seq.get(), {"core"}, {high});
        CHECK(element._get_ratios() == std::vector<double>{normal, high, normal});
    }

    SECTION("the weights track a body set that changed during setup") {
        auto seq = build("delete b1\n");
        REQUIRE(seq->_get_rigidbody()->molecule.size_body() == 2);

        RelativeHydrationElement element(seq.get(), {"b3"}, {high});
        CHECK(element._get_ratios() == std::vector<double>{normal, high});
    }

    SECTION("an unknown body name is rejected") {
        auto seq = build();
        CHECK_THROWS(RelativeHydrationElement(seq.get(), {"doesnotexist"}, {high}));
    }

    SECTION("a symmetry replica is rejected: it has no hydration of its own to scale") {
        auto seq = build("symmetry b1 c2\n");
        REQUIRE(seq->setup()._body_name_registry().contains("b1s1r1"));
        CHECK_THROWS(RelativeHydrationElement(seq.get(), {"b1s1r1"}, {high}));
    }

    SECTION("the script form parses and produces one weight per body") {
        auto seq = build("relative_hydration b2 max\n");
        REQUIRE(seq != nullptr);
        CHECK(seq->_get_rigidbody()->molecule.size_body() == 3);
    }
}
