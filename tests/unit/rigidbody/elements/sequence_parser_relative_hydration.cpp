// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>

#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/detail/BodyIndexOps.h>
#include <rigidbody/sequencer/elements/setup/RelativeHydrationElement.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/Rigidbody.h>
#include <data/Molecule.h>
#include <settings/All.h>
#include <io/ExistingFile.h>

#include <support/temp_file.h>

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
        test::TempFile config("ausaxs_seq_rel_hydration_test", ".conf", content);
        SequenceParser parser;
        return parser.parse_file(config);
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

        RelativeHydrationElement element(seq.get(), "b2", maximum);
        CHECK(element._get_ratios() == std::vector<double>{normal, maximum, normal});
    }

    SECTION("several named bodies each land at their own index") {
        auto seq = build();
        RelativeHydrationElement third(seq.get(), "b3", high);
        RelativeHydrationElement first(seq.get(), "b1", maximum);
        CHECK(first._get_ratios() == std::vector<double>{maximum, normal, high});
    }

    SECTION("naming every body leaves no default weights") {
        auto seq = build();
        RelativeHydrationElement e1(seq.get(), "b1", maximum);
        RelativeHydrationElement e2(seq.get(), "b2", high);
        RelativeHydrationElement e3(seq.get(), "b3", normal);
        CHECK(e3._get_ratios() == std::vector<double>{maximum, high, normal});
    }

    SECTION("a renamed body is addressable by its alias") {
        auto seq = build("rename b2 core\n");
        RelativeHydrationElement element(seq.get(), "core", high);
        CHECK(element._get_ratios() == std::vector<double>{normal, high, normal});
    }

    SECTION("a level declared through an alias survives a later rename") {
        auto seq = build("rename b2 core\n");
        RelativeHydrationElement element(seq.get(), "core", high);
        seq->setup()._body_name_registry().rename("core", "other");
        CHECK(element._get_ratios() == std::vector<double>{normal, high, normal});
    }

    SECTION("the weights track a body set that changed during setup") {
        auto seq = build("delete b1\n");
        REQUIRE(seq->_get_rigidbody()->molecule.size_body() == 2);

        RelativeHydrationElement element(seq.get(), "b3", high);
        CHECK(element._get_ratios() == std::vector<double>{normal, high});
    }

    SECTION("an unknown body name is rejected") {
        auto seq = build();
        CHECK_THROWS(RelativeHydrationElement(seq.get(), "doesnotexist", high));
    }

    SECTION("a symmetry replica is rejected: it has no hydration of its own to scale") {
        auto seq = build("symmetry b1 c2\n");
        REQUIRE(seq->setup()._body_name_registry().contains("b1s1r1"));
        CHECK_THROWS(RelativeHydrationElement(seq.get(), "b1s1r1", high));
    }

    SECTION("the script form parses and produces one weight per body") {
        auto seq = build("relative_hydration b2 max\n");
        REQUIRE(seq != nullptr);
        CHECK(seq->_get_rigidbody()->molecule.size_body() == 3);
    }
}

TEST_CASE_METHOD(SequenceParserRelativeHydrationFixture, "RelativeHydrationElement: declarations accumulate", "[files]") {
    // The ratio vector spans the whole body set, so an element cannot describe only its own bodies. Before the shared store,
    // a second element rebuilt the vector from scratch and reset every body the first one had named back to normal.
    constexpr double normal = 1.0, high = 1.5, maximum = 1.75, low = 0.5;

    SECTION("a second element does not reset the first one's bodies") {
        auto seq = build();
        RelativeHydrationElement first(seq.get(), "b1", maximum);
        RelativeHydrationElement second(seq.get(), "b3", low);
        CHECK(second._get_ratios() == std::vector<double>{maximum, normal, low});
        CHECK(first._get_ratios() == std::vector<double>{maximum, normal, low}); // both read the same store
    }

    SECTION("the last declaration for a body wins") {
        auto seq = build();
        RelativeHydrationElement first(seq.get(), "b2", maximum);
        RelativeHydrationElement second(seq.get(), "b2", low);
        CHECK(second._get_ratios() == std::vector<double>{normal, low, normal});
    }

    SECTION("naming one body two ways collapses onto a single entry") {
        auto seq = build("rename b2 core\n");
        RelativeHydrationElement byname(seq.get(), "b2", maximum);
        RelativeHydrationElement byalias(seq.get(), "core", low);
        CHECK(byalias._get_ratios() == std::vector<double>{normal, low, normal});
    }

    SECTION("the store does not leak into a later sequencer") {
        {
            auto seq = build();
            RelativeHydrationElement element(seq.get(), "b1", maximum);
            REQUIRE(element._get_ratios() == std::vector<double>{maximum, normal, normal});
        }
        auto seq = build();
        RelativeHydrationElement element(seq.get(), "b3", high);
        CHECK(element._get_ratios() == std::vector<double>{normal, normal, high});
    }

    SECTION("deleting a body that was given a level is an error, not a silent drop") {
        auto seq = build();
        RelativeHydrationElement element(seq.get(), "b1", maximum);
        sequencer::detail::erase_bodies(seq.get(), {0});
        CHECK_THROWS(element._get_ratios());
    }

    SECTION("repeated script lines accumulate, replacing the removed block form") {
        auto seq = build(
            "relative_hydration b1 max\n"
            "relative_hydration b3 low\n"
        );
        REQUIRE(seq != nullptr);

        // read the merged store back off one of the real elements the script produced
        RelativeHydrationElement* element = nullptr;
        for (const auto* list : {&seq->setup()._get_elements(), &seq->_get_elements()}) {
            for (const auto& e : *list) {
                if (auto* rh = dynamic_cast<RelativeHydrationElement*>(e.get())) {element = rh;}
            }
        }
        REQUIRE(element != nullptr);
        CHECK(element->_get_ratios() == std::vector<double>{maximum, normal, low});
    }
}
