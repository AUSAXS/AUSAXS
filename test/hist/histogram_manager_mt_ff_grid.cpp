#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hist/distance_calculator/HistogramManagerMT.h>
#include <hist/distance_calculator/HistogramManagerMTFFGrid.h>
#include <hist/distance_calculator/HistogramManagerMTFFGridSurface.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvg.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGrid.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGridSurface.h>
#include <data/Body.h>
#include <data/Molecule.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <grid/Grid.h>
#include <settings/All.h>
#include <utility/Utility.h>

#include "hist_test_helper.h"
#include "../grid/grid_debug.h"

using namespace hist;
using namespace data;
using namespace data::record;

auto test = [] (const Molecule& protein, std::function<std::unique_ptr<ICompositeDistanceHistogram>(const Molecule&)> calculate) {
    settings::molecule::center = false; // to avoid rounding errors
    auto h = calculate(protein);

    // convert the grid to water atoms
    auto exv_grid = protein.get_grid()->generate_excluded_volume(false);
    std::vector<Water> waters(exv_grid.interior.size());
    for (unsigned int i = 0; i < exv_grid.interior.size(); i++) {
        waters[i] = Water(i, "C", "", "LYS", 'A', 1, "", exv_grid.interior[i], 1, 0, constants::atom_t::dummy, "");
        waters[i].set_effective_charge(1);
    }

    Molecule exv(protein.get_atoms(), waters);
    auto h_exv = hist::HistogramManager<true>(&exv).calculate_all();

    // calculate the xx, ax, aa distributions
    auto h_cast = static_cast<CompositeDistanceHistogramFFAvg*>(h.get());
    hist::Distribution1D xx1, ax1, aa1;
    {
        // we have to sum all form factor contributions into the single distributions
        auto aa = h_cast->get_aa_counts_ff();
        hist::Distribution1D temp_aa(aa.size_z()), temp_ax(aa.size_z()), temp_xx(aa.size_z());
        for (unsigned int i = 0; i < form_factor::get_count_without_excluded_volume(); ++i) {
            for (unsigned int j = 0; j < form_factor::get_count_without_excluded_volume(); ++j) {
                std::transform(aa.begin(i, j), aa.end(i, j), temp_aa.begin(), temp_aa.begin(), std::plus<>());
            }
            std::transform(aa.begin(i, form_factor::exv_bin), aa.end(i, form_factor::exv_bin), temp_ax.begin(), temp_ax.begin(), std::plus<>());
            std::transform(aa.begin(form_factor::exv_bin, i), aa.end(form_factor::exv_bin, i), temp_ax.begin(), temp_ax.begin(), std::plus<>()); // should all be 0
        }
        std::transform(aa.begin(form_factor::exv_bin, form_factor::exv_bin), aa.end(form_factor::exv_bin, form_factor::exv_bin), temp_xx.begin(), temp_xx.begin(), std::plus<>());

        aa1 = std::move(temp_aa);
        ax1 = std::move(temp_ax);
        xx1 = std::move(temp_xx);
    }

    auto xx2 = h_exv->get_ww_counts();
    auto ax2 = h_exv->get_aw_counts();
    auto aa2 = h_exv->get_aa_counts();

    if (xx1.size() < xx2.size()) {
        for (unsigned int i = xx1.size(); i < xx2.size(); ++i) {
            REQUIRE(xx2.index(i) == 0);
            REQUIRE(ax2.index(i) == 0);
            REQUIRE(aa2.index(i) == 0);
        }
    } else {
        for (unsigned int i = xx2.size(); i < xx1.size(); ++i) {
            REQUIRE(xx1.index(i) == 0);
            REQUIRE(ax1.index(i) == 0);
            REQUIRE(aa1.index(i) == 0);
        }
    }

    // aa
    for (int i = 0; i < std::min<int>(xx1.size(), xx2.size()); ++i) {
        if (!utility::approx(aa1.index(i), aa2.index(i), 1e-6, 0)) {
            std::cout << "histogram_manager_mt_ff_grid failed at index " << i << std::endl;
            REQUIRE_THAT(aa1.index(i), Catch::Matchers::WithinAbs(aa2.index(i), 1e-6));
        }
    }
    SUCCEED();

    // xx
    for (int i = 0; i < std::min<int>(xx1.size(), xx2.size()); ++i) {
        if (!utility::approx(xx1.index(i), xx2.index(i), 1e-6, 0)) {
            std::cout << "histogram_manager_mt_ff_grid failed at index " << i << std::endl;
            REQUIRE_THAT(xx1.index(i), Catch::Matchers::WithinAbs(xx2.index(i), 1e-6));
        }
    } 
    SUCCEED();

    // ax
    for (int i = 0; i < std::min<int>(xx1.size(), xx2.size()); ++i) {
        if (!utility::approx(ax1.index(i), ax2.index(i), 1e-6, 0)) {
            std::cout << "histogram_manager_mt_ff_grid failed at index " << i << std::endl;
            REQUIRE_THAT(ax1.index(i), Catch::Matchers::WithinAbs(ax2.index(i), 1e-6));
        }
    }
    SUCCEED();
};

// Check that the Grid histograms are correct
TEST_CASE("HistogramManagerMTFFGrid::calculate", "[files]") {
    settings::molecule::use_effective_charge = false;
    SECTION("simple") {
        settings::grid::width = GENERATE(0.2, 0.5, 1, 2);
        settings::grid::exv_radius = settings::grid::width;

        SECTION(std::string("width = ") + std::to_string(settings::grid::width)) {
            Atom a1(0, "C", "", "LYS", 'A', 1, "", {0, 0, 0}, 1, 0, constants::atom_t::dummy, "");
            Molecule protein({a1});
            test(protein, [](const Molecule& protein) {return hist::HistogramManagerMTFFGrid(&protein).calculate_all();});
        }
    }

    SECTION("actual data") {
        settings::general::verbose = false;
        settings::grid::width = 1;
        settings::grid::exv_radius = 1;
        Molecule protein("test/files/LAR1-2.pdb");
        protein.clear_hydration();
        test(protein, [](const Molecule& protein) {return hist::HistogramManagerMTFFGrid(&protein).calculate_all();});
    }
}

// Check that the GridSurface histograms are correct
TEST_CASE("HistogramManagerMTFFGridSurface::calculate", "[files]") {
    settings::molecule::use_effective_charge = false;
    settings::molecule::center = false;
    settings::general::verbose = false;

    SECTION("simple") {
        settings::grid::width = GENERATE(0.2, 0.5, 1, 2);
        settings::grid::exv_radius = settings::grid::width;

        SECTION(std::string("width = ") + std::to_string(settings::grid::width)) {
            Atom a1(0, "C", "", "LYS", 'A', 1, "", {0, 0, 0}, 1, 0, constants::atom_t::dummy, "");
            Molecule protein({a1});
            test(protein, [](const Molecule& protein) {return hist::HistogramManagerMTFFGridSurface(&protein).calculate_all();});
        }
    }

    SECTION("actual data") {
        settings::general::verbose = false;
        settings::grid::width = 1;
        settings::grid::exv_radius = 1;
        Molecule protein("test/files/LAR1-2.pdb");
        protein.clear_hydration();
        test(protein, [](const Molecule& protein) {return hist::HistogramManagerMTFFGridSurface(&protein).calculate_all();});
    }

    SECTION("compare with Grid") {
        settings::grid::min_bins = 20;

        Molecule protein("test/files/LAR1-2.pdb");
        protein.generate_new_hydration();

        auto h_grid  = hist::HistogramManagerMTFFGrid(&protein).calculate_all();
        auto h_grids = hist::HistogramManagerMTFFGridSurface(&protein).calculate_all();

        auto h_grid_cast = static_cast<CompositeDistanceHistogramFFGrid*>(h_grid.get());
        auto aa1 = h_grid_cast->get_aa_counts_ff();
        auto ax1 = h_grid_cast->get_aw_counts_ff();
        auto xx1 = h_grid_cast->get_ww_counts_ff();

        auto h_grids_cast = static_cast<CompositeDistanceHistogramFFGridSurface*>(h_grids.get());
        auto aa2 = h_grids_cast->get_aa_counts_ff();
        auto ax2 = h_grids_cast->get_aw_counts_ff();
        auto xx2 = h_grids_cast->get_ww_counts_ff();

        CHECK(xx1.size() == xx2.size());
        for (unsigned int k = 0; k < xx1.size(); ++k) {
            if (!utility::approx(xx1.index(k), xx2.index(k), 1e-3, 0)) {
                std::cout << "histogram_manager_mt_ff_grid failed at index " << k << std::endl;
                REQUIRE_THAT(xx2.index(k), Catch::Matchers::WithinAbs(xx1.index(k), 1e-3));
            }
            SUCCEED();
        }

        CHECK((ax1.size_x() == ax2.size_x() && ax1.size_y() == ax2.size_y()));
        for (unsigned int k = 0; k < ax2.size_y(); ++k) {
            for (unsigned int j = 0; j < ax2.size_x(); ++j) {
                if (!utility::approx(ax2.index(j, k), ax1.index(j, k), 1e-3, 0)) {
                    std::cout << "histogram_manager_mt_ff_grid failed at index " << j << ", " << k << std::endl;
                    REQUIRE_THAT(ax2.index(j, k), Catch::Matchers::WithinAbs(ax1.index(j, k), 1e-3));
                }
                SUCCEED();
            }
        }

        CHECK((aa1.size_x() == aa2.size_x() && aa1.size_y() == aa2.size_y() && aa1.size_z() == aa2.size_z()));
        for (unsigned int k = 0; k < aa2.size_z(); ++k) {
            for (unsigned int j = 0; j < aa2.size_y(); ++j) {
                for (unsigned int i = 0; i < aa2.size_x(); ++i) {
                    if (!utility::approx(aa2.index(i, j, k), aa1.index(i, j, k), 1e-3, 0)) {
                        std::cout << "histogram_manager_mt_ff_grid failed at index " << i << ", " << j << ", " << k << std::endl;
                        REQUIRE_THAT(aa2.index(i, j, k), Catch::Matchers::WithinAbs(aa1.index(i, j, k), 1e-3));
                    }
                    SUCCEED();
                }
            }
        }
    }
}

// Check that the weighted bins are correct and separate for the excluded volume and the protein atoms
TEST_CASE("HistogramManagerMTFFGrid: weighted_bins", "[files]") {
    settings::molecule::use_effective_charge = false;
    settings::hist::weighted_bins = true;
    Molecule protein("test/files/2epe.pdb");

    SECTION("simple") {
        std::vector<Atom> b1 = {Atom(Vector3<double>(-1, -1, -1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b2 = {Atom(Vector3<double>( 1, -1, -1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>( 1, 1, -1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b3 = {Atom(Vector3<double>(-1, -1,  1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>(-1, 1,  1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Atom> b4 = {Atom(Vector3<double>( 1, -1,  1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>( 1, 1,  1), 1, constants::atom_t::C, "C", 1)};
        std::vector<Body> a = {Body(b1), Body(b2), Body(b3), Body(b4)};
        Molecule protein(a);
        set_unity_charge(protein);

        auto exv_grid = protein.get_grid()->generate_excluded_volume(false);
        std::vector<Atom> atoms(exv_grid.interior.size());
        for (unsigned int i = 0; i < exv_grid.interior.size(); i++) {
            atoms[i] = Atom(i, "C", "", "LYS", 'A', 1, "", exv_grid.interior[i], 1, 0, constants::atom_t::dummy, "");
        }
        Molecule exv(atoms);

        auto h_grid  = hist::HistogramManagerMTFFGrid(&protein).calculate_all();
        auto h_grids = hist::HistogramManagerMTFFGridSurface(&protein).calculate_all();
        auto h_exv   = hist::HistogramManagerMT<true>(&exv).calculate_all();
        auto h_atom  = hist::HistogramManagerMT<true>(&protein).calculate_all();

        auto h_grid_cast = static_cast<CompositeDistanceHistogramFFGrid*>(h_grid.get());
        auto h_grids_cast = static_cast<CompositeDistanceHistogramFFGridSurface*>(h_grids.get());

        CHECK(compare_hist(h_grid_cast->get_d_axis(),   h_atom->get_d_axis()));
        CHECK(compare_hist(h_grids_cast->get_d_axis(),  h_atom->get_d_axis()));

        CHECK(compare_hist(h_grid_cast->get_d_axis_xx(), h_exv->get_d_axis()));
        CHECK(compare_hist(h_grids_cast->get_d_axis_xx(), h_exv->get_d_axis()));
    }

    SECTION("simple, all") {
        std::vector<Atom> atoms = SimpleCube::atoms;
        atoms.push_back(Atom(Vector3<double>(0, 0, 0), 1, constants::atom_t::C, "C", 1));
        std::for_each(atoms.begin(), atoms.end(), [](Atom& a) {a.set_effective_charge(1);});

        Molecule protein(atoms);
        GridDebug::generate_debug_grid(protein); // overrides exv generation
        auto h = hist::HistogramManagerMTFFGrid(&protein).calculate_all();
        auto h_cast = static_cast<CompositeDistanceHistogramFFGrid*>(h.get());

        // check the distance axes
        REQUIRE(SimpleCube::check_exact(h_cast->get_d_axis()));
        REQUIRE(SimpleCube::check_exact(h_cast->get_d_axis_ax()));
        REQUIRE(SimpleCube::check_exact(h_cast->get_d_axis_xx()));
    }

    SECTION("real data") {
        auto exv_grid = protein.get_grid()->generate_excluded_volume(false);
        std::vector<Atom> atoms(exv_grid.interior.size());
        for (unsigned int i = 0; i < exv_grid.interior.size(); i++) {
            atoms[i] = Atom(i, "C", "", "LYS", 'A', 1, "", exv_grid.interior[i], 1, 0, constants::atom_t::dummy, "");
        }
        Molecule exv(atoms);
        for (auto& b : exv.get_bodies()) {for (auto& a : b.get_atoms()) {a.set_effective_charge(1);}}

        auto h_grid  = hist::HistogramManagerMTFFGrid(&protein).calculate_all();
        auto h_grids = hist::HistogramManagerMTFFGridSurface(&protein).calculate_all();
        auto h_exv   = hist::HistogramManagerMT<true>(&exv).calculate_all();
        auto h_atom  = hist::HistogramManagerMT<true>(&protein).calculate_all();

        auto h_grid_cast = static_cast<CompositeDistanceHistogramFFGrid*>(h_grid.get());
        auto h_grids_cast = static_cast<CompositeDistanceHistogramFFGridSurface*>(h_grid.get());

        CHECK(compare_hist(h_grid_cast->get_d_axis(),   h_atom->get_d_axis()));
        CHECK(compare_hist(h_grids_cast->get_d_axis(),   h_atom->get_d_axis()));

        CHECK(compare_hist(h_grid_cast->get_d_axis_xx(), h_exv->get_d_axis()));
        CHECK(compare_hist(h_grids_cast->get_d_axis_xx(), h_exv->get_d_axis()));
    }
}

auto calc_scat = [] (double k) {
    const auto& q_axis = constants::axes::q_vals;
    auto ff_C = form_factor::storage::atomic::get_form_factor(form_factor::form_factor_t::C);

    auto V = std::pow(2*settings::grid::exv_radius, 3);
    form_factor::FormFactor ffx = form_factor::ExvFormFactor(V);
    auto d = SimpleCube::d_exact;

    std::vector<double> Iq_exp(q_axis.size(), 0);
    for (unsigned int q = 0; q < q_axis.size(); ++q) {
        // calculation: 8 points (scaled by k)
        //          1 line  of length 0
        //          3 lines of length 2
        //          3 lines of length sqrt(2*2^2) = sqrt(8) = 2.82
        //          1 line  of length sqrt(3*2^2) = sqrt(12) = 3.46
        //
        // calculation: 1 center point
        //          1 line  of length 0
        //          16 lines of length sqrt(3) = 1.73 (counting both directions)

        double aasum = 
            9 + 
            16*std::sin(q_axis[q]*d[1])/(q_axis[q]*d[1]) + 
            24*std::sin(q_axis[q]*d[2])/(q_axis[q]*d[2]) + 
            24*std::sin(q_axis[q]*d[3])/(q_axis[q]*d[3]) + 
            8*std::sin(q_axis[q]*d[4])/(q_axis[q]*d[4]);
        double axsum = 
            1 +
            k*8 +
            8*(1+k)*std::sin(q_axis[q]*d[1])/(q_axis[q]*d[1]) + 
            k*24*std::sin(q_axis[q]*d[2])/(q_axis[q]*d[2]) + 
            k*24*std::sin(q_axis[q]*d[3])/(q_axis[q]*d[3]) + 
            k*8*std::sin(q_axis[q]*d[4])/(q_axis[q]*d[4]);
        double xxsum = 
            1 +
            k*k*8 + 
            k*16*std::sin(q_axis[q]*d[1])/(q_axis[q]*d[1]) + 
            k*k*24*std::sin(q_axis[q]*d[2])/(q_axis[q]*d[2]) + 
            k*k*24*std::sin(q_axis[q]*d[3])/(q_axis[q]*d[3]) + 
            k*k*8*std::sin(q_axis[q]*d[4])/(q_axis[q]*d[4]);

        // if (q==0) {
        //     std::cout << "aasum: " << aasum << std::endl;
        //     std::cout << "\t9*1" << std::endl;
        //     std::cout << "\t16*sin(" << q_axis[q] << "*" << d[1] << ")/(" << q_axis[q] << "*" << d[1] << ") = " << std::sin(q_axis[q]*d[1])/(q_axis[q]*d[1]) << std::endl;
        //     std::cout << "\t24*sin(" << q_axis[q] << "*" << d[2] << ")/(" << q_axis[q] << "*" << d[2] << ") = " << std::sin(q_axis[q]*d[2])/(q_axis[q]*d[2]) << std::endl;
        //     std::cout << "\t24*sin(" << q_axis[q] << "*" << d[3] << ")/(" << q_axis[q] << "*" << d[3] << ") = " << std::sin(q_axis[q]*d[3])/(q_axis[q]*d[3]) << std::endl;
        //     std::cout << "\t8*sin(" << q_axis[q] << "*" << d[4] << ")/(" << q_axis[q] << "*" << d[4] << ") = " << std::sin(q_axis[q]*d[4])/(q_axis[q]*d[4]) << std::endl;

        //     std::cout << "axsum: " << axsum << std::endl;
        //     std::cout << "\t1" << std::endl;
        //     std::cout << "\t" << 8*k << std::endl;
        //     std::cout << "\t" << 24*k << "*sin(" << q_axis[q] << "*" << d[2] << ")/(" << q_axis[q] << "*" << d[2] << ") = " << std::sin(q_axis[q]*d[2])/(q_axis[q]*d[2]) << std::endl;
        //     std::cout << "\t" << 24*k << "*sin(" << q_axis[q] << "*" << d[3] << ")/(" << q_axis[q] << "*" << d[3] << ") = " << std::sin(q_axis[q]*d[3])/(q_axis[q]*d[3]) << std::endl;
        //     std::cout << "\t" << 8*k << " *sin(" << q_axis[q] << "*" << d[4] << ")/(" << q_axis[q] << "*" << d[4] << ") = " << std::sin(q_axis[q]*d[4])/(q_axis[q]*d[4]) << std::endl;

        //     std::cout << "xxsum: " << xxsum << std::endl;
        //     std::cout << "\t1" << std::endl;
        //     std::cout << "\t" << 8*k*k << std::endl;
        //     std::cout << "\t" << 16 << "*sin(" << q_axis[q] << "*" << d[1] << ")/(" << q_axis[q] << "*" << d[1] << ") = " << k*k*std::sin(q_axis[q]*d[1])/(q_axis[q]*d[1]) << std::endl;
        //     std::cout << "\t" << 24*k*k << "*sin(" << q_axis[q] << "*" << d[2] << ")/(" << q_axis[q] << "*" << d[2] << ") = " << std::sin(q_axis[q]*d[2])/(q_axis[q]*d[2]) << std::endl;
        //     std::cout << "\t" << 24*k*k << "*sin(" << q_axis[q] << "*" << d[3] << ")/(" << q_axis[q] << "*" << d[3] << ") = " << std::sin(q_axis[q]*d[3])/(q_axis[q]*d[3]) << std::endl;
        //     std::cout << "\t" << 8*k*k << "*sin(" << q_axis[q] << "*" << d[4] << ")/(" << q_axis[q] << "*" << d[4] << ") = " << std::sin(q_axis[q]*d[4])/(q_axis[q]*d[4]) << std::endl;
        // }

        Iq_exp[q] += aasum*std::pow(ff_C.evaluate(q_axis[q]), 2);               // + aa
        Iq_exp[q] -= 2*axsum*ff_C.evaluate(q_axis[q])*ffx.evaluate(q_axis[q]);  // -2ax
        Iq_exp[q] += xxsum*std::pow(ffx.evaluate(q_axis[q]), 2);                // + xx
    }
    return Iq_exp;
};

// Check that the surface scaling works as expected
TEST_CASE("HistogramManagerMTFFGridSurface: surface_scaling") {
    settings::molecule::use_effective_charge = false;
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = true;

    settings::grid::save_exv = true;
    settings::general::output = "temp/test/hist/hmmtffg/";
    settings::grid::rvol = 2;
    std::vector<Atom> atoms = SimpleCube::atoms;
    atoms.push_back(Atom(Vector3<double>(0, 0, 0), 1, constants::atom_t::C, "C", 1));
    std::for_each(atoms.begin(), atoms.end(), [](Atom& a) {a.set_effective_charge(1);});

    Molecule protein(atoms);
    GridDebug::generate_debug_grid(protein);
    auto h = hist::HistogramManagerMTFFGridSurface(&protein).calculate_all();
    auto h_cast = static_cast<CompositeDistanceHistogramFFGridSurface*>(h.get());

    // check the distance axes
    REQUIRE(SimpleCube::check_exact(h_cast->get_d_axis()));
    REQUIRE(SimpleCube::check_exact(h_cast->get_d_axis_ax()));
    REQUIRE(SimpleCube::check_exact(h_cast->get_d_axis_xx()));

    // x1
    auto Iq_exp = calc_scat(1);
    REQUIRE(compare_hist(Iq_exp, h->debye_transform()));

    // x2
    Iq_exp = calc_scat(2);
    h_cast->apply_excluded_volume_scaling_factor(2);
    REQUIRE(compare_hist(Iq_exp, h->debye_transform()));

    // x3
    Iq_exp = calc_scat(3);
    h_cast->apply_excluded_volume_scaling_factor(3);
    REQUIRE(compare_hist(Iq_exp, h->debye_transform()));

    // x0.5
    Iq_exp = calc_scat(0.5);
    h_cast->apply_excluded_volume_scaling_factor(0.5);
    REQUIRE(compare_hist(Iq_exp, h->debye_transform()));
}