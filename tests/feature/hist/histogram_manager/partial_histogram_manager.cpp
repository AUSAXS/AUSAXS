#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/state/Signaller.h>  // IWYU pragma: keep
#include <hist/histogram_manager/HistogramManagerMT.h>
#include <hist/histogram_manager/PartialHistogramManager.h>
#include <hist/histogram_manager/PartialHistogramManagerMT.h>
#include <hist/histogram_manager/PartialSymmetryManagerMT.h>
#include <settings/All.h>

#include <numeric>
#include <random>

#include <support/hist_test_helper.h>

using namespace ausaxs;
using namespace ausaxs::data;

// Test that the first calculation is correct
TEST_CASE("PartialHistogramManager: initial calculation", "[gpu]") {
    settings::exv::exv_method = settings::exv::ExvMethod::None; // the expected histograms use the unmodified atomic weights
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::grid::min_bins = 100;
    std::vector<std::string> files = {
        "2epe",
        "6lyz",
        "c60",
        "diamond",
        "LAR1-2"
    };
    for (const auto& f : files) {
        {   // no hydration
            data::Molecule protein("tests/files/" + f + ".pdb");
            protein.clear_hydration();
            auto p_exp = hist::HistogramManagerMT<true>(&protein).calculate_all()->get_weighted_counts();
            {   // phm
                auto phm = hist::PartialHistogramManager<true>(&protein).calculate_all()->get_weighted_counts();
                REQUIRE(compare_hist_approx(p_exp, phm, 0, 1e-2));
            }
            {   // phm_mt
                auto phm_mt = hist::PartialHistogramManagerMT<true>(&protein).calculate_all()->get_weighted_counts();
                REQUIRE(compare_hist_approx(p_exp, phm_mt, 0, 1e-2));
            }
            {   // pshm_mt
                auto pshm_mt = hist::PartialSymmetryManagerMT<true>(&protein).calculate_all()->get_weighted_counts();
                REQUIRE(compare_hist_approx(p_exp, pshm_mt, 0, 1e-2));
            }
        }

        {   // with hydration
            data::Molecule protein("tests/files/" + f + ".pdb");
            protein.generate_new_hydration();
            auto p_exp = hist::HistogramManagerMT<true>(&protein).calculate_all()->get_weighted_counts();
            {   // phm
                auto phm = hist::PartialHistogramManager<true>(&protein).calculate_all()->get_weighted_counts();
                REQUIRE(compare_hist_approx(p_exp, phm, 0, 1e-2));
            }
            {   // phm_mt
                auto phm_mt = hist::PartialHistogramManagerMT<true>(&protein).calculate_all()->get_weighted_counts();
                REQUIRE(compare_hist_approx(p_exp, phm_mt, 0, 1e-2));
            }
            {   // pshm_mt
                auto pshm_mt = hist::PartialSymmetryManagerMT<true>(&protein).calculate_all()->get_weighted_counts();
                REQUIRE(compare_hist_approx(p_exp, pshm_mt, 0, 1e-2));
            }
        }
    }
}

static auto test = [] (data::Molecule& protein, auto&& phm) {
    // no changes
    auto p_exp = hist::HistogramManagerMT<true>(&protein).calculate_all()->get_weighted_counts();
    auto phm_res = phm(protein)->get_weighted_counts();
    REQUIRE(compare_hist_approx(p_exp, phm_res, 0, 1e-2));

    // change hydration
    protein.generate_new_hydration();
    p_exp = hist::HistogramManagerMT<true>(&protein).calculate_all()->get_weighted_counts();
    phm_res = phm(protein)->get_weighted_counts();
    REQUIRE(compare_hist_approx(p_exp, phm_res, 0, 1e-2));

    // external change
    protein.get_body(1).translate({1, 1, 1});
    p_exp = hist::HistogramManagerMT<true>(&protein).calculate_all()->get_weighted_counts();
    phm_res = phm(protein)->get_weighted_counts();
    REQUIRE(compare_hist_approx(p_exp, phm_res, 0, 1e-2));

    // internal change
    protein.get_body(1).get_atom(0).weight() = 100;
    protein.get_body(1).get_signaller()->modified_internal();
    p_exp = hist::HistogramManagerMT<true>(&protein).calculate_all()->get_weighted_counts();
    phm_res = phm(protein)->get_weighted_counts();
    REQUIRE(compare_hist_approx(p_exp, phm_res, 0, 1e-2));
};

static auto test_random = [] (data::Molecule& protein, auto&& phm) {
    auto p_exp = hist::HistogramManagerMT<true>(&protein).calculate_all()->get_weighted_counts();
    auto phm_res = phm(protein)->get_weighted_counts();
    REQUIRE(compare_hist_approx(p_exp, phm_res, 0, 1e-2));

    bool modify_external = GENERATE(false, true);
    bool modify_internal = GENERATE(false, true);
    bool modify_hydration = GENERATE(false, true);

    static std::random_device seed;
    static std::mt19937 gen(seed());
    static std::uniform_int_distribution<> ri(0, 100);
    static std::uniform_real_distribution<> rd(-10, 10);

    for (int i = 0; i < 2; ++i) {
        if (modify_external) {
            int body_index = ri(gen) % protein.size_body();
            protein.get_body(body_index).translate({rd(gen), rd(gen), rd(gen)});
        }

        if (modify_internal) {
            int body_index = ri(gen) % protein.size_body();
            int atom_index = ri(gen) % protein.get_body(body_index).size_atom();
            for (int j = 0; j < ri(gen); ++j) {
                protein.get_body(body_index).get_atom(atom_index).weight() = std::abs(rd(gen));
            }
            protein.get_body(body_index).get_signaller()->modified_internal();
        }

        if (modify_hydration) {
            protein.generate_new_hydration();
        }

        auto p_exp = hist::HistogramManagerMT<true>(&protein).calculate_all()->get_weighted_counts();
        auto phm_res = phm(protein)->get_weighted_counts();
        REQUIRE(compare_hist_approx(p_exp, phm_res, 0, 1e-2));
    }
};

// Test that subsequent calculations are correct
TEST_CASE("PartialHistogramManager: subsequent calculations", "[gpu]") {
    settings::exv::exv_method = settings::exv::ExvMethod::None; // the expected histograms use the unmodified atomic weights
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    data::Molecule protein({
        Body("tests/files/2epe.pdb"), 
        Body{std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)}}
    });

    protein.generate_new_hydration();

    // deterministic
    test(protein, [] (const Molecule& protein) {return hist::PartialHistogramManager<true>(&protein).calculate_all();});
    test(protein, [] (const Molecule& protein) {return hist::PartialHistogramManagerMT<true>(&protein).calculate_all();});
    test(protein, [] (const Molecule& protein) {return hist::PartialSymmetryManagerMT<true>(&protein).calculate_all();});    

    // random
    test_random(protein, [] (const Molecule& protein) {return hist::PartialHistogramManager<true>(&protein).calculate_all();});
    test_random(protein, [] (const Molecule& protein) {return hist::PartialHistogramManagerMT<true>(&protein).calculate_all();});
    test_random(protein, [] (const Molecule& protein) {return hist::PartialSymmetryManagerMT<true>(&protein).calculate_all();});    
}

TEST_CASE("PartialHistogramManager: grows its axis when the structure outgrows it", "[gpu]") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    auto choice = GENERATE(
        settings::hist::HistogramManagerChoice::PartialHistogramManager,
        settings::hist::HistogramManagerChoice::PartialHistogramManagerMT,
        settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT
    );

    data::Molecule protein({Body("tests/files/2epe.pdb"), Body("tests/files/6lyz.pdb")});
    protein.generate_new_hydration();
    protein.set_histogram_manager(choice);
    auto before = protein.get_histogram()->get_weighted_counts();

    protein.get_body(1).translate({300, 0, 0});
    auto after = protein.get_histogram()->get_weighted_counts();

    REQUIRE_THAT(std::reduce(after.begin(), after.end(), 0.0), Catch::Matchers::WithinRel(std::reduce(before.begin(), before.end(), 0.0), 1e-9));
    REQUIRE(after.size() > before.size());
}

TEST_CASE("PartialHistogramManager: a total-only calculation does not leave stale partials behind", "[gpu]") {
    settings::exv::exv_method = settings::exv::ExvMethod::None; // the expected histograms use the unmodified atomic weights
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    auto choice = GENERATE(
        settings::hist::HistogramManagerChoice::PartialHistogramManager,
        settings::hist::HistogramManagerChoice::PartialHistogramManagerMT
    );

    data::Molecule protein({Body("tests/files/2epe.pdb"), Body("tests/files/6lyz.pdb")});
    protein.generate_new_hydration();
    protein.set_histogram_manager(choice);
    (void) protein.get_histogram();

    // a total-only calculation of a modified state, as OverlapConstraint::evaluate() does, followed by a full one
    protein.get_body(1).translate({20, 0, 0});
    protein.generate_new_hydration();
    (void) protein.get_total_histogram();
    auto h = protein.get_histogram();

    auto ref = hist::HistogramManagerMT<true>(&protein).calculate_all();
    REQUIRE(h->get_aa_counts().size() == h->get_weighted_counts().size());
    CHECK(compare_hist_approx(h->get_aa_counts(), ref->get_aa_counts()));
    CHECK(compare_hist_approx(h->get_aw_counts(), ref->get_aw_counts()));
    CHECK(compare_hist_approx(h->get_ww_counts(), ref->get_ww_counts()));
}

TEST_CASE("PartialHistogramManager: the simple excluded volume follows the current geometry", "[gpu]") {
    settings::general::verbose = false;
    settings::exv::exv_method = settings::exv::ExvMethod::Simple;
    auto choice = GENERATE(
        settings::hist::HistogramManagerChoice::PartialHistogramManager,
        settings::hist::HistogramManagerChoice::PartialHistogramManagerMT,
        settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT
    );

    // the first two bodies overlap, so pulling one away changes the grid volume and with it the displaced charge per atom
    Molecule protein({Body("tests/files/2epe.pdb"), Body("tests/files/2epe.pdb"), Body("tests/files/6lyz.pdb")});
    protein.get_body(1).translate({5, 0, 0});
    protein.get_body(2).translate({60, 0, 0});
    protein.set_histogram_manager(choice);
    auto I = [&] () {return protein.get_histogram()->debye_transform().get_intensity();};
    auto I0 = I();
    double V0 = protein.get_volume_grid();

    protein.get_body(0).translate({30, 0, 0});
    protein.clear_grid();
    (void) I();
    REQUIRE(protein.get_volume_grid() != V0);

    // body 2 is rebuilt while the volume differs. Once everything is back in place, the result must not depend on that
    protein.get_body(2).translate({5, 0, 0});
    (void) I();
    protein.get_body(2).translate({-5, 0, 0});
    (void) I();
    protein.get_body(0).translate({-30, 0, 0});
    protein.clear_grid();
    REQUIRE(protein.get_volume_grid() == V0);

    auto I1 = I();
    REQUIRE(I1.size() == I0.size());
    for (unsigned int i = 0; i < I0.size(); ++i) {
        CHECK_THAT(I1[i], Catch::Matchers::WithinRel(I0[i], 1e-5));
    }
}
