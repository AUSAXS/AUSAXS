#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hist/distance_calculator/HistogramManagerMTFFGrid.h>
#include <hist/distance_calculator/HistogramManagerMTFFGridSurface.h>
#include <hist/distance_calculator/HistogramManagerMTFFGridScalableExv.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGrid.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGridSurface.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGridScalableExv.h>
#include <data/Molecule.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <grid/Grid.h>
#include <dataset/SimpleDataset.h>
#include <plots/All.h>
#include <settings/All.h>

#include "hist/hist_test_helper.h"
#include "grid/grid_debug.h"

using namespace data;
using namespace data::record;

TEST_CASE("CompositeDistanceHistogramFFGrid::volumes", "[manual]") {
    settings::general::verbose = false;
    data::Molecule protein("tests/files/2epe.pdb");
    auto V = protein.get_volume_grid();

    std::vector<double> volumes;
    std::vector<double> rxs;
    for (double rx = 0.1; rx < 1; rx += 0.1) {
        settings::grid::cell_width = rx;
        settings::grid::exv::width = rx;
        protein.clear_grid();
        volumes.push_back(protein.get_volume_grid());
        rxs.push_back(rx);

        // hist::CompositeDistanceHistogramFFGrid::regenerate_table();
        // auto h = hist::HistogramManagerMTFFGrid<false>(&protein).calculate_all();
        // auto h_cast = static_cast<hist::CompositeDistanceHistogramFFGrid*>(h.get());
        // auto ffx = h_cast->get_ff_table().index(form_factor::exv_bin, form_factor::exv_bin);

        // auto gridx = protein.get_grid()->generate_excluded_volume();
        // volumes.push_back(gridx.size()*ffx.evaluate(0)/constants::charge::density::water);
        // rxs.push_back(rx);

        // CHECK_THAT(V, Catch::Matchers::WithinRel(protein.get_volume_grid(), 1e-1));
    }

    SimpleDataset dataset(rxs, volumes);
    plots::PlotDataset::quick_plot(dataset, plots::PlotOptions({{"xlabel", "Grid width [Å]"}, {"ylabel", "Volume [Å³]"}, {"color", style::color::blue}}), "temp/tests/hist/composite_distance_histogram_ff_grid/volumes.png");
}

TEST_CASE("CompositeDistanceHistogramFFGridSurface: compare_profiles") {}


auto calc_scat = [] (double k) {
    const auto& q_axis = constants::axes::q_vals;
    auto ff_C = form_factor::storage::atomic::get_form_factor(form_factor::form_factor_t::C);

    auto V = std::pow(2*settings::grid::exv::width, 3);
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

        auto cx = hist::CompositeDistanceHistogramFFGridSurface::exv_factor(constants::axes::q_vals[q], k);
        double aasum = 
            9 + 
            16*std::sin(q_axis[q]*d[1])/(q_axis[q]*d[1]) + 
            24*std::sin(q_axis[q]*d[2])/(q_axis[q]*d[2]) + 
            24*std::sin(q_axis[q]*d[3])/(q_axis[q]*d[3]) + 
            8*std::sin(q_axis[q]*d[4])/(q_axis[q]*d[4]);
        double axsum = 
            1 +
            cx*8 +
            8*(1+cx)*std::sin(q_axis[q]*d[1])/(q_axis[q]*d[1]) + 
            cx*24*std::sin(q_axis[q]*d[2])/(q_axis[q]*d[2]) + 
            cx*24*std::sin(q_axis[q]*d[3])/(q_axis[q]*d[3]) + 
            cx*8*std::sin(q_axis[q]*d[4])/(q_axis[q]*d[4]);
        double xxsum = 
            1 +
            cx*cx*8 + 
            cx*16*std::sin(q_axis[q]*d[1])/(q_axis[q]*d[1]) + 
            cx*cx*24*std::sin(q_axis[q]*d[2])/(q_axis[q]*d[2]) + 
            cx*cx*24*std::sin(q_axis[q]*d[3])/(q_axis[q]*d[3]) + 
            cx*cx*8*std::sin(q_axis[q]*d[4])/(q_axis[q]*d[4]);

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

auto calc_scat_water = [] () {
    const auto& q_axis = constants::axes::q_vals;
    auto ff_C = form_factor::storage::atomic::get_form_factor(form_factor::form_factor_t::C);
    auto ff_O = form_factor::storage::atomic::get_form_factor(static_cast<form_factor::form_factor_t>(form_factor::water_bin));

    auto V = std::pow(2*settings::grid::exv::width, 3);
    form_factor::FormFactor ffx = form_factor::ExvFormFactor(V);
    auto d = SimpleCube::d_exact;

    std::vector<double> Iq_exp(q_axis.size(), 0);
    for (unsigned int q = 0; q < q_axis.size(); ++q) {
        // due to the simple setup, all distances are the same
        double aasum = 
            9 + 
            16*std::sin(q_axis[q]*d[1])/(q_axis[q]*d[1]) + 
            24*std::sin(q_axis[q]*d[2])/(q_axis[q]*d[2]) + 
            24*std::sin(q_axis[q]*d[3])/(q_axis[q]*d[3]) + 
            8*std::sin(q_axis[q]*d[4])/(q_axis[q]*d[4]);
        Iq_exp[q] += aasum*std::pow(ff_C.evaluate(q_axis[q]), 2);               // + aa
        Iq_exp[q] -= 2*aasum*ff_C.evaluate(q_axis[q])*ffx.evaluate(q_axis[q]);  // -2ax
        Iq_exp[q] += aasum*std::pow(ffx.evaluate(q_axis[q]), 2);                // + xx
        Iq_exp[q] += 2*aasum*ff_O.evaluate(q_axis[q])*ff_C.evaluate(q_axis[q]); // +2aw
        Iq_exp[q] -= 2*aasum*ff_O.evaluate(q_axis[q])*ffx.evaluate(q_axis[q]);  // -2wx
        Iq_exp[q] += aasum*std::pow(ff_O.evaluate(q_axis[q]), 2);               // + ww
    }
    return Iq_exp;
};

// Check that the debye transform is correct
TEST_CASE("HistogramManagerMTFFGrid::debye_transform") {
    settings::molecule::use_effective_charge = false;
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = true;
    settings::hist::weighted_bins = true;
    settings::grid::min_exv_radius = 1;
    settings::grid::exv::width = 0.5;
    settings::grid::min_exv_radius = 0;

    std::vector<Atom> atoms = SimpleCube::get_atoms();
    atoms.push_back(Atom(Vector3<double>(0, 0, 0), 1, constants::atom_t::C, "C", 1));
    std::for_each(atoms.begin(), atoms.end(), [](Atom& a) {a.set_effective_charge(1);});

    std::vector<Water> waters = SimpleCube::get_waters();
    waters.push_back(Water::create_new_water(Vector3<double>(0, 0, 0)));
    std::for_each(waters.begin(), waters.end(), [](Water& w) {w.set_effective_charge(1);});

    Molecule protein(atoms);
    GridDebug::generate_debug_grid(protein);
    protein.get_waters() = waters;

    SECTION("Grid") {
        auto h = hist::HistogramManagerMTFFGrid(&protein).calculate_all();
        auto h_cast = static_cast<hist::CompositeDistanceHistogramFFGrid*>(h.get());
        REQUIRE(SimpleCube::check_exact(h_cast->get_d_axis()));
        REQUIRE(SimpleCube::check_exact(h_cast->get_d_axis_ax()));
        REQUIRE(SimpleCube::check_exact(h_cast->get_d_axis_xx()));

        auto Iq_exp = calc_scat_water();
        REQUIRE(compare_hist(Iq_exp, h->debye_transform()));
    }

    SECTION("GridSurface") {
        auto h = hist::HistogramManagerMTFFGridSurface(&protein).calculate_all();
        auto h_cast = static_cast<hist::CompositeDistanceHistogramFFGridSurface*>(h.get());
        REQUIRE(SimpleCube::check_exact(h_cast->get_d_axis()));
        REQUIRE(SimpleCube::check_exact(h_cast->get_d_axis_ax()));
        REQUIRE(SimpleCube::check_exact(h_cast->get_d_axis_xx()));

        auto Iq_exp = calc_scat_water();
        REQUIRE(compare_hist(Iq_exp, h->debye_transform()));
    }

    SECTION("GridScalableExv") {
        auto h = hist::HistogramManagerMTFFGridScalableExv(&protein).calculate_all();
        auto h_cast = static_cast<hist::CompositeDistanceHistogramFFGridScalableExv*>(h.get());
        REQUIRE(SimpleCube::check_exact(h_cast->get_d_axis()));
        REQUIRE(SimpleCube::check_exact(h_cast->get_d_axis_ax()));
        REQUIRE(SimpleCube::check_exact(h_cast->get_d_axis_xx()));

        auto Iq_exp = calc_scat_water();
        REQUIRE(compare_hist(Iq_exp, h->debye_transform()));
    }
}

// Check that the surface scaling works as expected
TEST_CASE("HistogramManagerMTFFGridSurface: surface_scaling") {
    settings::molecule::use_effective_charge = false;
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = true;
    settings::hist::weighted_bins = true;
    settings::grid::cell_width = 1;
    settings::grid::exv::width = 0.5;
    settings::grid::min_exv_radius = 0;

    std::vector<Atom> atoms = SimpleCube::get_atoms();
    atoms.push_back(Atom(Vector3<double>(0, 0, 0), 1, constants::atom_t::C, "C", 1));
    std::for_each(atoms.begin(), atoms.end(), [](Atom& a) {a.set_effective_charge(1);});

    Molecule protein(atoms);
    GridDebug::generate_debug_grid(protein);
    auto h = hist::HistogramManagerMTFFGridSurface(&protein).calculate_all();
    auto h_cast = static_cast<hist::CompositeDistanceHistogramFFGridSurface*>(h.get());

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

    // x0
    Iq_exp = calc_scat(0);
    h_cast->apply_excluded_volume_scaling_factor(0);
    REQUIRE(compare_hist(Iq_exp, h->debye_transform()));
}