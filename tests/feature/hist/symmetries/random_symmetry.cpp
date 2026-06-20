#include <catch2/catch_test_macros.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <data/symmetry/CyclicSymmetry.h>
#include <data/symmetry/PointSymmetry.h>
#include <data/symmetry/TetrahedralSymmetry.h>
#include <data/symmetry/OctahedralSymmetry.h>
#include <data/symmetry/IcosahedralSymmetry.h>
#include <data/symmetry/CompositeSymmetry.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <hist/histogram_manager/SymmetryManagerMT.h>
#include <hist/histogram_manager/PartialSymmetryManagerMT.h>
#include <settings/All.h>

#include "hist/hist_test_helper.h"
#include "settings/HistogramSettings.h"

#include <numbers>
#include <random>

using namespace ausaxs;
using namespace ausaxs::data;

// Fuzz test: stack random symmetries of every type onto small random structures and verify the reuse-based SymmetryManagers 
// reproduce the brute-force explicit_structure histogram. The generators respect the calculator's hard limit on the 
// per-partial-histogram scale factor: at most x30 on any job, except exactly x60 for a standalone icosahedral group. The 
// dominant scale is the self-correlation, 1 + size_symmetry_total() per body, so each body's total copy count is bounded below:
//   - cyclic/point/composite stacks keep the per-body copy sum well under 29 (self-scale <= 30)
//   - polyhedral groups are applied standalone: T/O/I give self-scales 12/24/60 (60 is the special case)
namespace {
    std::mt19937& rng() {static std::mt19937 g(std::random_device{}()); return g;}
    double rnd(double lo, double hi) {return std::uniform_real_distribution<>(lo, hi)(rng());}
    int rndi(int lo, int hi) {return std::uniform_int_distribution<>(lo, hi)(rng());}
    Vector3<double> rnd_vec(double m) {return {rnd(-m, m), rnd(-m, m), rnd(-m, m)};}
    Vector3<double> rnd_axis() {
        Vector3<double> v = rnd_vec(1);
        if (v.magnitude() < 1e-3) {return {0, 0, 1};}
        v.normalize();
        return v;
    }

    // a cyclic c_n with n in [2, max_n]; contributes (n-1) copies, internal scale <= n
    std::unique_ptr<symmetry::ISymmetry> rnd_cyclic(int max_n) {
        int n = rndi(2, max_n);
        return std::make_unique<symmetry::CyclicSymmetry>(
            symmetry::CyclicSymmetry::_Relation{rnd_vec(8)},
            symmetry::CyclicSymmetry::_Repeat{rnd_axis(), 2*std::numbers::pi/n},
            n-1
        );
    }

    std::unique_ptr<symmetry::ISymmetry> rnd_point() {
        return std::make_unique<symmetry::PointSymmetry>(rnd_vec(8), rnd_vec(std::numbers::pi));
    }

    // a 2-level composite of small cyclics; product of orders <= 9, so self-scale and internal scales <= 9
    std::unique_ptr<symmetry::ISymmetry> rnd_composite() {
        return std::make_unique<symmetry::CompositeSymmetry>(rnd_cyclic(3), rnd_cyclic(3));
    }

    // a randomly-oriented polyhedral group; icosahedral only when it may stand alone (self-scale 60)
    std::unique_ptr<symmetry::ISymmetry> rnd_polyhedral(bool allow_icosahedral) {
        std::unique_ptr<symmetry::IPolyhedralSymmetry> p;
        switch (rndi(0, allow_icosahedral ? 2 : 1)) {
            case 0:  p = std::make_unique<symmetry::TetrahedralSymmetry>(); break;
            case 1:  p = std::make_unique<symmetry::OctahedralSymmetry>();  break;
            default: p = std::make_unique<symmetry::IcosahedralSymmetry>(); break;
        }
        p->translation = rnd_vec(8);
        p->rotation = rnd_vec(std::numbers::pi);
        return p;
    }

    // give one body a random symmetry configuration that stays within the scale limits
    void add_random_symmetries(Body& body) {
        switch (rndi(0, 4)) {
            case 0: body.symmetry().add(rnd_cyclic(6)); break;                 // single cyclic, self-scale <= 6
            case 1: body.symmetry().add(rnd_point()); break;                   // single point, self-scale 2
            case 2: body.symmetry().add(rnd_polyhedral(true)); break;          // standalone T/O/I, self-scale 12/24/60
            case 3: body.symmetry().add(rnd_composite()); break;               // single composite, self-scale <= 9
            case 4: {                                                          // a stack of light cyclic/point symmetries
                int budget = 10;                                               // keep self-scale = 1 + sum(copies) <= 11
                for (int k = 0; k < 3; ++k) {
                    auto s = rndi(0, 1) ? rnd_cyclic(3) : rnd_point();
                    int copies = static_cast<int>(s->repetitions());
                    if (copies > budget) {break;}
                    budget -= copies;
                    body.symmetry().add(std::move(s));
                }
                if (body.size_symmetry() == 0) {body.symmetry().add(rnd_point());}
                break;
            }
        }
    }

    void run_fuzz(settings::hist::HistogramManagerChoice choice) {
        for (int iter = 0; iter < 30; ++iter) {
            int n_bodies = rndi(1, 2);
            std::vector<Body> bodies;
            for (int b = 0; b < n_bodies; ++b) {
                std::vector<AtomFF> atoms;
                for (int j = 0, na = rndi(1, 3); j < na; ++j) {atoms.push_back(AtomFF(rnd_vec(8), form_factor::form_factor_t::C));}
                bodies.emplace_back(std::move(atoms));
            }
            Molecule m(std::move(bodies));
            m.set_histogram_manager(choice);
            for (int b = 0; b < n_bodies; ++b) {add_random_symmetries(m.get_body(b));}
            set_unity_charge(m);

            auto h = m.get_histogram()->get_weighted_counts();

            // ground truth: materialise every copy explicitly and histogram with the plain manager
            std::vector<Body> explicit_bodies;
            for (int b = 0; b < n_bodies; ++b) {
                auto e = m.get_body(b).symmetry().explicit_structure();
                explicit_bodies.emplace_back(std::move(e.atoms), std::move(e.waters));
            }
            Molecule m2(std::move(explicit_bodies));
            set_unity_charge(m2);
            auto h2 = m2.get_histogram()->get_weighted_counts();

            // equalise lengths so compare_hist_approx's +-1 window is not truncated at the tail: a reuse representative and the 
            // individually-binned explicit pairs of the same group can differ by sub-bin FP (their relative transforms match only 
            // to the bucketer tolerance), and near a bin edge at the histogram tail that benign 1-bin drift would otherwise escape
            // the window where one histogram has already been trimmed to its last non-zero bin.
            std::size_t n = std::max(h.size(), h2.size());
            h.resize(n, 0); h2.resize(n, 0);

            INFO("fuzz iteration " << iter);
            CHECK(compare_hist_approx(h, h2));
        }
    }
}

TEST_CASE("SymmetryManager: random mixed-symmetry fuzz") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    SECTION("SymmetryManager") {
        run_fuzz(settings::hist::HistogramManagerChoice::HistogramSymmetryManagerMT);
    }
    SECTION("PartialSymmetryManager") {
        run_fuzz(settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT);
    }
}
