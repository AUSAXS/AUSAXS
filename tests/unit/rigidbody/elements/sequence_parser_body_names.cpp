// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>

#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/detail/BodyNameRegistry.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/elements/setup/BodySymmetrySelector.h>
#include <rigidbody/Rigidbody.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <settings/All.h>
#include <io/ExistingFile.h>

#include <support/temp_file.h>

#include <algorithm>
#include <set>
#include <string>
#include <vector>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::sequencer;
namespace seqdetail = ausaxs::rigidbody::sequencer::detail; // "detail" alone is ambiguous between the sequencer, rigidbody, and ausaxs namespaces

struct SequenceParserBodyNamesFixture {
    SequenceParserBodyNamesFixture() {
        settings::general::verbose = false;
        settings::molecule::implicit_hydrogens = false;
        settings::grid::min_bins = 100;
        settings::hydrate::hydration_strategy = settings::hydrate::HydrationStrategy::NoStrategy;
        settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    }

    std::unique_ptr<Sequencer> parse(const std::string& content) {
        test::TempFile config(".conf", content);
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

    // Every registered entity (base body or symmetry replica) must be reachable under its own name: no two of them may share a default name, or one of the
    // two is left addressable by nobody.
    static void expect_unique_names(const seqdetail::BodyNameRegistry& registry) {
        std::set<std::string> seen;
        for (const auto& [index, entry] : registry.all()) {
            INFO("duplicate default name \"" << entry.default_name << "\"");
            CHECK(seen.insert(entry.default_name).second);
        }
    }

    // Each name in `names` must resolve to a base body, and no two of them to the same one.
    static void expect_distinct_bodies(const seqdetail::BodyNameRegistry& registry, const std::vector<std::string>& names) {
        std::set<int> seen;
        for (const auto& name : names) {
            INFO("body name \"" << name << "\"");
            REQUIRE(registry.contains(name));
            CHECK(seen.insert(registry.resolve_body(name)).second);
        }
    }

    // Lowest residue sequence id carried by the body `name` refers to, used to tell the fragments of a split apart by the residues they hold.
    static int first_residue(Sequencer& seq, const std::string& name) {
        const auto& body = seq._get_rigidbody()->molecule.get_body(seq.setup()._get_body(name));
        const auto& metadata = body.get_metadata();
        REQUIRE(metadata.has_value());
        REQUIRE(metadata->residue_seq.has_value());
        return *std::min_element(metadata->residue_seq->begin(), metadata->residue_seq->end());
    }
};

TEST_CASE_METHOD(SequenceParserBodyNamesFixture, "BodyNameRegistry: default names stay unique as bodies are created and destroyed", "[files]") {
    SECTION("a plain load yields exactly b1..bN") {
        auto seq = build();
        REQUIRE(seq != nullptr);
        REQUIRE(seq->_get_rigidbody()->molecule.size_body() == 3);
        CHECK(seq->setup()._body_name_registry().base_body_names() == std::vector<std::string>{"b1", "b2", "b3"});
    }

    SECTION("splitting a body mid-script does not reuse a surviving body's name") {
        // regression: the fragments used to be named from their new indices, regenerating "b3" - a name the third body already held and kept through the
        // reindexing, leaving the fragment at residues 40-59 unaddressable
        auto seq = build("split b2 60\n");
        REQUIRE(seq != nullptr);
        auto& registry = seq->setup()._body_name_registry();

        REQUIRE(seq->_get_rigidbody()->molecule.size_body() == 4);
        auto names = registry.base_body_names();
        REQUIRE(names.size() == 4);
        expect_unique_names(registry);
        expect_distinct_bodies(registry, names);

        // the leading fragment continues b2's identity, so "b2" survives the split; only the trailing fragment is a new body. Note the fragments are
        // appended at the tail, so they come after b3 in index order.
        CHECK(names == std::vector<std::string>{"b1", "b3", "b2", "b4"});
        CHECK(first_residue(*seq, "b2") == 40);
        CHECK(first_residue(*seq, "b4") == 60);
    }

    SECTION("a split body's alias follows its leading fragment") {
        auto seq = build("rename b2 core\nsplit core 60\n");
        REQUIRE(seq != nullptr);
        auto& registry = seq->setup()._body_name_registry();

        REQUIRE(seq->_get_rigidbody()->molecule.size_body() == 4);
        expect_unique_names(registry);
        expect_distinct_bodies(registry, registry.base_body_names());

        // both of the original's names still address the same body: the fragment holding the residues the original started with
        REQUIRE(registry.contains("core"));
        CHECK(registry.resolve_body("core") == registry.resolve_body("b2"));
        CHECK(first_residue(*seq, "core") == 40);
        CHECK(registry.base_body_names() == std::vector<std::string>{"b1", "b3", "core", "b4"});
    }

    SECTION("a body created after a delete collides with nothing, and the survivors keep their names") {
        auto seq = build("delete b1\nsplit b3 100\n");
        REQUIRE(seq != nullptr);
        auto& registry = seq->setup()._body_name_registry();

        REQUIRE(seq->_get_rigidbody()->molecule.size_body() == 3);
        auto names = registry.base_body_names();
        REQUIRE(names.size() == 3);
        expect_unique_names(registry);
        expect_distinct_bodies(registry, names);

        // b2 survived the delete and the split untouched, b3 followed its leading fragment, and the trailing fragment continues the sequence past b3
        CHECK_FALSE(registry.contains("b1"));
        CHECK(names == std::vector<std::string>{"b2", "b3", "b4"});
        CHECK(first_residue(*seq, "b3") == 80);
        CHECK(first_residue(*seq, "b4") == 100);
    }

    SECTION("a deleted body's name is never reissued") {
        // a script referring to "b2" after a delete-then-create must fail loudly rather than silently retarget the new body
        auto seq = build("delete b2\ncopy b1 c1\ncopy b3 c2\n");
        REQUIRE(seq != nullptr);
        auto& registry = seq->setup()._body_name_registry();

        REQUIRE(seq->_get_rigidbody()->molecule.size_body() == 4);
        expect_unique_names(registry);
        CHECK_FALSE(registry.contains("b2"));
        expect_distinct_bodies(registry, registry.base_body_names());
    }

    SECTION("a fragment never picks up a deleted body's name") {
        // the trailing fragment draws from the counter, which must still skip past the retired b2 rather than reissuing it
        auto seq = build("delete b2\nsplit b1 20\n");
        REQUIRE(seq != nullptr);
        auto& registry = seq->setup()._body_name_registry();

        REQUIRE(seq->_get_rigidbody()->molecule.size_body() == 3);
        expect_unique_names(registry);
        CHECK_FALSE(registry.contains("b2"));
        CHECK(registry.base_body_names() == std::vector<std::string>{"b3", "b1", "b4"});
    }

    SECTION("a copied body continues the default-name sequence") {
        auto seq = build("copy b1 clone\n");
        REQUIRE(seq != nullptr);
        auto& registry = seq->setup()._body_name_registry();

        REQUIRE(seq->_get_rigidbody()->molecule.size_body() == 4);
        expect_unique_names(registry);
        // the copy is addressable both by its alias and by its own default name, which must not shadow any of b1..b3
        CHECK(registry.base_body_names() == std::vector<std::string>{"b1", "b2", "b3", "clone"});
        CHECK(registry.contains("b4"));
        CHECK(registry.resolve_body("b4") == registry.resolve_body("clone"));
    }
}

TEST_CASE_METHOD(SequenceParserBodyNamesFixture, "BodyNameRegistry: replica tags follow their own body's default name", "[files]") {
    // b2 carries a symmetry and is then split, so the fragments are created after b3 already exists: their replica tags must be built from the fragments' own
    // (post-split) default names rather than from their indices, which no longer identify them
    auto seq = build("symmetry b2 c2\nsplit b2 60\n");
    REQUIRE(seq != nullptr);
    auto& registry = seq->setup()._body_name_registry();
    REQUIRE(seq->_get_rigidbody()->molecule.size_body() == 4);

    expect_unique_names(registry);
    expect_distinct_bodies(registry, registry.base_body_names());

    // the leading fragment inherited b2's name, so the tags of its replicas are unchanged by the split too
    CHECK(registry.contains("b2s1r1"));
    CHECK(registry.resolve_body("b2") == registry.resolve("b2s1r1").body);

    int replicas = 0;
    for (const auto& [index, entry] : registry.all()) {
        auto sel = seqdetail::from_index(index);
        if (sel.replica == 0) {continue;} // base body, not a replica tag
        ++replicas;

        // every replica tag is prefixed by the default name of the body it belongs to
        const std::string& base = registry.entry(seqdetail::to_index(sel.body)).default_name;
        INFO("replica tag \"" << entry.default_name << "\" of body \"" << base << "\"");
        CHECK(entry.default_name.rfind(base + "s", 0) == 0);
    }
    CHECK(0 < replicas);
}
