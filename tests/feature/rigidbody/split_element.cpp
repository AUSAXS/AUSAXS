// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/elements/All.h>
#include <rigidbody/Rigidbody.h>
#include <data/symmetry/ReferenceSymmetry.h>
#include <data/symmetry/CompositeSymmetry.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <settings/All.h>
#include <io/ExistingFile.h>

#include <hist/hist_test_helper.h>

#include <algorithm>
#include <fstream>
#include <random>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::sequencer;

namespace {
    std::unique_ptr<Sequencer> parse(const std::string& content) {
        static int counter = 0;
        std::string path = "/tmp/ausaxs_split_element_test_" + std::to_string(counter++) + ".conf";
        std::ofstream f(path);
        f << content;
        f.close();
        SequenceParser parser;
        return parser.parse_file(path);
    }

    // Load 2epe as a single body, optionally applying `symmetry_name` to it (empty = none), then run `extra` script
    // lines (e.g. a split). Returns the fully-built rigidbody. 2epe is a small, real, multi-element structure: the
    // mass-weighted combined centre of mass used by ReferenceSymmetry only reproduces the pre-split single-body
    // centre of mass exactly if elements of different mass are actually split across fragments, so a uniform-mass
    // toy structure would not exercise this.
    std::unique_ptr<Sequencer> build(const std::string& symmetry_name, const std::string& extra = "") {
        std::string script =
            "load {\n"
            "    pdb tests/files/2epe.pdb\n"
            "    saxs tests/files/2epe.dat\n"
            "}\n";
        if (!symmetry_name.empty()) {script += "symmetry b1 " + symmetry_name + "\n";}
        return parse(script + extra);
    }

    // n_splits equidistant residue ids spanning `body`'s own residue range, so the caller does not depend on
    // hardcoded knowledge of the structure's numbering.
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

    // Translate every body of `molecule` by the identical vector. Moving a whole symmetric assembly this way must
    // leave its scattering unchanged whether it is represented as one body carrying its own symmetry, or as several
    // fragments sharing a ReferenceSymmetry - which is exactly what SplitElement is supposed to preserve.
    //
    // Deliberately translation-only, not rotation: a body's own orientation and its symmetry's axis/frame are
    // independent, separately-optimizable parameters in this design (confirmed against ReferenceSymmetry's own
    // ground-truth test) - rotating a body's atoms directly does not rotate the fixed axis its symmetry stores, so
    // an arbitrary external rotation is not expected to reproduce the same assembly either way, split or not. This
    // is a pre-existing property of the symmetry model, unrelated to Split; translation alone already exercises the
    // exact invariant Split depends on (ReferenceSymmetry::combined_cm tracking the original body's centre of mass).
    void rigid_translate(data::Molecule& molecule, const Vector3<double>& t) {
        for (std::size_t i = 0; i < molecule.size_body(); ++i) {molecule.get_body(i).translate(t);}
    }

    // Assert that `unsplit` and `split` scatter identically, both immediately and after a handful of random rigid
    // translations of the whole assembly.
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
    // isolate the comparison to atom positions and the shared symmetry; explicit hydration is placed by a
    // randomized, grid-dependent procedure that would otherwise differ between the two independently-hydrated
    // molecules regardless of whether Split is correct
    settings::hydrate::hydration_strategy = settings::hydrate::HydrationStrategy::NoStrategy;

    // cover the full range of bases the shared ReferenceSymmetry now accepts: cyclic (c2/c3), dihedral (d2, a
    // single polyhedral leaf) and a genuine composite (p2-c3, two leaves) - the composite in particular exercises
    // the for_each_leaf descent through the ReferenceSymmetry wrapper that this feature adds
    auto symmetry_name = GENERATE(as<std::string>{}, "c2", "c3", "d2", "p2-c3");
    int expected_leaves = std::string(symmetry_name).find('-') == std::string::npos ? 1 : 2;
    int n_splits = GENERATE(1, 2, 3);

    auto seq_unsplit = build(symmetry_name); // kept alive: it owns the rigidbody below
    auto rb_unsplit = seq_unsplit->_get_rigidbody();
    REQUIRE(rb_unsplit != nullptr);
    REQUIRE(rb_unsplit->molecule.size_body() == 1);

    // 2epe.pdb carries crystallographic waters, loaded as explicit hydration on b1. Split has no way to assign an
    // existing water molecule to one fragment over another, so its fragments are always built waterless (see
    // partition_by_residue) - clear the unsplit reference's waters too, so both sides start from the same
    // (water-free) state and the comparison isolates the atom+symmetry logic Split is actually responsible for.
    rb_unsplit->molecule.get_body(0).clear_hydration();

    auto ids = equidistant_split_ids(rb_unsplit->molecule.get_body(0), n_splits);
    auto seq_split = build(symmetry_name, split_line(ids));
    auto rb_split = seq_split->_get_rigidbody();
    REQUIRE(rb_split != nullptr);
    REQUIRE(rb_split->molecule.size_body() == static_cast<std::size_t>(n_splits+1));

    // Split must tie every fragment into one shared symmetric assembly: the first fragment owns a ReferenceSymmetry,
    // the rest hold non-owning views onto it. (A regression that instead gave each fragment an independent copy of
    // the symmetry would already break the scattering comparison below, but assert the intended topology directly.)
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
