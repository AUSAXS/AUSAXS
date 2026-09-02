// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/elements/All.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/constraints/IDistanceConstraint.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <rigidbody/detail/MoleculeTransformParametersAbsolute.h>
#include <rigidbody/parameters/BodyTransformParametersAbsolute.h>
#include <rigidbody/parameters/ParameterGenerationStrategy.h>
#include <rigidbody/transform/TransformStrategy.h>
#include <data/symmetry/ReferenceSymmetry.h>
#include <data/symmetry/CompositeSymmetry.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <settings/All.h>
#include <io/ExistingFile.h>

#include <hist/hist_test_helper.h>

#include <algorithm>
#include <random>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::sequencer;

namespace {
    std::unique_ptr<Sequencer> parse(const std::string& content) {
        return SequenceParser().parse_text(content);
    }

    // Load 2epe as a single body, optionally applying `symmetry_name` to it (empty = none), then run `extra` script lines (e.g. a split).
    std::unique_ptr<Sequencer> build(const std::string& symmetry_name, const std::string& extra = "") {
        std::string script =
            "load {\n"
            "    pdb tests/files/2epe.pdb\n"
            "    saxs tests/files/2epe.dat\n"
            "}\n";
        if (!symmetry_name.empty()) {script += "symmetry b1 " + symmetry_name + "\n";}
        return parse(script + extra);
    }

    // n_splits equidistant residue ids spanning `body`'s own residue range, so the caller does not depend on hardcoded knowledge of the structure's numbering.
    std::vector<int> equidistant_split_ids(const data::Body& body, int n_splits) {
        const auto& metadata = body.get_metadata();
        REQUIRE(metadata.has_value());
        REQUIRE(metadata->residue_seq.has_value());
        const auto& resseq = *metadata->residue_seq;
        int min_id = *std::min_element(resseq.begin(), resseq.end());
        int max_id = *std::max_element(resseq.begin(), resseq.end());

        std::vector<int> ids;
        for (int k = 1; k <= n_splits; ++k) {ids.push_back(min_id + (max_id - min_id)*k/(n_splits+1));}
        return ids;
    }

    std::string split_line(const std::vector<int>& ids) {
        std::string line = "split b1";
        for (int id : ids) {line += " " + std::to_string(id);}
        return line + "\n";
    }

    // Translate every body of `molecule` by the identical vector. Moving a whole symmetric assembly this way must leave its scattering unchanged whether it is 
    // represented as one body carrying its own symmetry, or as several fragments sharing a ReferenceSymmetry. 
    void rigid_translate(data::Molecule& molecule, const Vector3<double>& t) {
        for (std::size_t i = 0; i < molecule.size_body(); ++i) {molecule.get_body(i).translate(t);}
    }

    // Assert that `unsplit` and `split` scatter identically, both immediately and after a handful of random rigid translations of the whole assembly.
    void expect_same_scattering(data::Molecule& unsplit, data::Molecule& split) {
        CHECK(compare_hist_approx(unsplit.get_histogram()->get_weighted_counts(), split.get_histogram()->get_weighted_counts()));

        static std::random_device rd;
        static std::mt19937 gen(rd());
        std::uniform_real_distribution<> d_t(-5, 5);
        for (int i = 0; i < 5; ++i) {
            Vector3<double> t{d_t(gen), d_t(gen), d_t(gen)};
            rigid_translate(unsplit, t);
            rigid_translate(split, t);
            CHECK(compare_hist_approx(unsplit.get_histogram()->get_weighted_counts(), split.get_histogram()->get_weighted_counts()));
        }
    }
}

TEST_CASE("SplitElement: splitting a symmetric body preserves scattering under rigid transforms", "[files]") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::grid::min_bins = 100;

    // isolate the comparison to atom positions and the shared symmetry
    settings::hydrate::hydration_strategy = settings::hydrate::HydrationStrategy::NoStrategy;
    auto symmetry_name = GENERATE(as<std::string>{}, "c2", "c3", "d2", "p2-c3");
    int expected_leaves = std::string(symmetry_name).find('-') == std::string::npos ? 1 : 2;
    int n_splits = GENERATE(1, 2, 3);

    auto seq_unsplit = build(symmetry_name); // kept alive: it owns the rigidbody below
    auto rb_unsplit = seq_unsplit->_get_rigidbody();
    REQUIRE(rb_unsplit != nullptr);
    REQUIRE(rb_unsplit->molecule.size_body() == 1);
    rb_unsplit->molecule.get_body(0).clear_hydration();

    auto ids = equidistant_split_ids(rb_unsplit->molecule.get_body(0), n_splits);
    auto seq_split = build(symmetry_name, split_line(ids));
    auto rb_split = seq_split->_get_rigidbody();
    REQUIRE(rb_split != nullptr);
    REQUIRE(rb_split->molecule.size_body() == static_cast<std::size_t>(n_splits+1));

    // Split must tie every fragment into one shared symmetric assembly: the first fragment owns a ReferenceSymmetry, the rest hold non-owning views onto it. 
    auto* ref = dynamic_cast<symmetry::ReferenceSymmetry*>(rb_split->molecule.get_body(0).symmetry().get(0));
    REQUIRE(ref != nullptr);
    for (std::size_t b = 1; b < rb_split->molecule.size_body(); ++b) {
        auto* view = dynamic_cast<symmetry::ReferenceSymmetryView*>(rb_split->molecule.get_body(b).symmetry().get(0));
        REQUIRE(view != nullptr);
        CHECK(view->target() == ref);
    }

    // for_each_leaf must see through the ReferenceSymmetry into its base's own leaves rather than stopping at the wrapper
    std::vector<symmetry::ISymmetry*> leaves;
    symmetry::for_each_leaf(*ref, [&](symmetry::ISymmetry& leaf) {leaves.push_back(&leaf);});
    CHECK(static_cast<int>(leaves.size()) == expected_leaves);

    expect_same_scattering(rb_unsplit->molecule, rb_split->molecule);
}

TEST_CASE("SplitElement: splitting a body with no symmetry yields independent fragments", "[files]") {
    settings::flags::max_bin_count = constants::axes::d_axis.bins;
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::grid::min_bins = 100;
    settings::hydrate::hydration_strategy = settings::hydrate::HydrationStrategy::NoStrategy;

    int n_splits = GENERATE(1, 2, 3);

    auto seq_unsplit = build(""); // kept alive: it owns the rigidbody below
    auto rb_unsplit = seq_unsplit->_get_rigidbody();
    REQUIRE(rb_unsplit != nullptr);
    REQUIRE(rb_unsplit->molecule.size_body() == 1);
    rb_unsplit->molecule.get_body(0).clear_hydration();
    // the split molecule is (re)built on a PartialSymmetryManagerMT; match it on the unsplit reference so the
    // comparison isolates the atom partitioning rather than any difference between histogram-manager backends
    rb_unsplit->molecule.set_histogram_manager(settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT);

    auto ids = equidistant_split_ids(rb_unsplit->molecule.get_body(0), n_splits);
    auto seq_split = build("", split_line(ids));
    auto rb_split = seq_split->_get_rigidbody();
    REQUIRE(rb_split != nullptr);
    REQUIRE(rb_split->molecule.size_body() == static_cast<std::size_t>(n_splits+1));

    // no shared symmetry is created when the base body carries none
    for (std::size_t b = 0; b < rb_split->molecule.size_body(); ++b) {
        CHECK(rb_split->molecule.get_body(b).size_symmetry() == 0);
    }

    // partitioning the atoms into separate bodies must not change the total scattering
    expect_same_scattering(rb_unsplit->molecule, rb_split->molecule);
}

TEST_CASE("SplitElement: constrained optimization steps of split symmetric fragments", "[files]") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::grid::min_bins = 100;
    settings::hydrate::hydration_strategy = settings::hydrate::HydrationStrategy::NoStrategy;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Backbone;
    settings::rigidbody::transform_strategy = GENERATE(
        settings::rigidbody::TransformationStrategyChoice::RigidTransform,
        settings::rigidbody::TransformationStrategyChoice::SingleTransform
    );

    auto seq_probe = build("c2"); // only needed to read the residue numbering off the unsplit body
    auto ids = equidistant_split_ids(seq_probe->_get_rigidbody()->molecule.get_body(0), 3);

    auto seq = build("c2", split_line(ids) + "autoconstrain backbone\n");
    auto rb = seq->_get_rigidbody();
    REQUIRE(rb != nullptr);
    REQUIRE(rb->molecule.size_body() == 4);

    // the fragments do not all carry the same kind of symmetry - the primary owns a ReferenceSymmetry while the others hold ReferenceSymmetryViews - so a delta
    // generated for one of them can only be applied to that very body, and never to whichever body the transformed branch happens to begin with
    for (unsigned int ibody = 0; ibody < rb->molecule.size_body(); ++ibody) {
        const auto& constraints = rb->constraints->get_body_constraints(ibody);
        REQUIRE(!constraints.empty());

        auto shared_parameters = [&] {
            // re-resolved every time: undo() move-assigns into the body, replacing its symmetry storage
            auto* ref = dynamic_cast<symmetry::ReferenceSymmetry*>(rb->molecule.get_body(0).symmetry().get(0));
            REQUIRE(ref != nullptr);
            auto t = ref->span_translation();
            auto r = ref->span_rotation();
            std::vector<double> values(t.begin(), t.end());
            values.insert(values.end(), r.begin(), r.end());
            return values;
        };
        auto before = shared_parameters();

        auto par = rb->parameter_generator->next(ibody);
        REQUIRE(par.symmetry_pars.has_value());
        REQUIRE(par.symmetry_pars->size() == rb->molecule.get_body(ibody).size_symmetry());
        rb->transformer->apply(std::move(par), constraints[0], ibody);

        // the shared symmetry is driven through its owner alone; a delta generated for one of the views must leave it untouched
        auto after = shared_parameters();
        REQUIRE(after.size() == before.size());
        bool changed = false;
        for (std::size_t i = 0; i < after.size(); ++i) {changed |= after[i] != before[i];}
        CHECK(changed == (ibody == 0));

        rb->transformer->undo();
        auto restored = shared_parameters();
        REQUIRE(restored.size() == before.size());
        for (std::size_t i = 0; i < restored.size(); ++i) {CHECK(restored[i] == before[i]);}
    }
}
