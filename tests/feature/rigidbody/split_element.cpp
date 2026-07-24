// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/elements/All.h>
#include <rigidbody/Rigidbody.h>
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
}

TEST_CASE("SplitElement: splitting a symmetric body preserves scattering under rigid transforms", "[files]") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::grid::min_bins = 100;
    // isolate the comparison to atom positions and the shared symmetry; explicit hydration is placed by a
    // randomized, grid-dependent procedure that would otherwise differ between the two independently-hydrated
    // molecules regardless of whether Split is correct
    settings::hydrate::hydration_strategy = settings::hydrate::HydrationStrategy::NoStrategy;

    // a small, real, multi-element structure - the mass-weighted combined centre of mass used by ReferenceSymmetry
    // only reproduces the pre-split single-body centre of mass exactly if elements of different mass are actually
    // split across fragments, so a uniform-mass toy structure would not exercise this
    auto symmetry_name = GENERATE(as<std::string>{}, "c2", "c3", "d2");
    int n_splits = GENERATE(1, 2, 3);

    std::string script =
        "load {\n"
        "    pdb tests/files/2epe.pdb\n"
        "    saxs tests/files/2epe.dat\n"
        "}\n"
        "symmetry b1 " + symmetry_name + "\n";

    auto seq_unsplit = parse(script);
    auto rb_unsplit = seq_unsplit->_get_rigidbody();
    REQUIRE(rb_unsplit != nullptr);
    REQUIRE(rb_unsplit->molecule.size_body() == 1);

    // 2epe.pdb carries crystallographic waters, loaded as explicit hydration on b1. Split has no way to assign an
    // existing water molecule to one fragment over another, so its fragments are always built waterless (see
    // partition_by_residue) - clear the unsplit reference's waters too, so both sides start from the same
    // (water-free) state and the comparison isolates the atom+symmetry logic Split is actually responsible for.
    rb_unsplit->molecule.get_body(0).clear_hydration();

    // pick n_splits equidistant residue ids spanning the loaded body's own residue range, so the test does not
    // depend on hardcoded knowledge of the structure's numbering
    const auto& metadata = rb_unsplit->molecule.get_body(0).get_metadata();
    REQUIRE(metadata.has_value());
    REQUIRE(metadata->residue_seq.has_value());
    const auto& resseq = *metadata->residue_seq;
    int min_id = *std::min_element(resseq.begin(), resseq.end());
    int max_id = *std::max_element(resseq.begin(), resseq.end());

    std::string split_line = "split b1";
    for (int k = 1; k <= n_splits; ++k) {
        split_line += " " + std::to_string(min_id + (max_id - min_id)*k/(n_splits+1));
    }

    auto seq_split = parse(script + split_line + "\n");
    auto rb_split = seq_split->_get_rigidbody();
    REQUIRE(rb_split != nullptr);
    REQUIRE(rb_split->molecule.size_body() == static_cast<std::size_t>(n_splits+1));

    auto h0_unsplit = rb_unsplit->molecule.get_histogram()->get_weighted_counts();
    auto h0_split   = rb_split->molecule.get_histogram()->get_weighted_counts();
    CHECK(compare_hist_approx(h0_unsplit, h0_split));

    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<> d_t(-5, 5);

    for (int i = 0; i < 5; ++i) {
        Vector3<double> t{d_t(gen), d_t(gen), d_t(gen)};

        rigid_translate(rb_unsplit->molecule, t);
        rigid_translate(rb_split->molecule, t);

        auto h_unsplit = rb_unsplit->molecule.get_histogram()->get_weighted_counts();
        auto h_split   = rb_split->molecule.get_histogram()->get_weighted_counts();
        CHECK(compare_hist_approx(h_unsplit, h_split));
    }
}
