#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hist/histogram_manager/HistogramManagerMTFFGrid.h>
#include <hist/histogram_manager/HistogramManagerMTFFGridSurface.h>
#include <hist/histogram_manager/HistogramManagerMTFFGridScalableExv.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGrid.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGridSurface.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGridScalableExv.h>
#include <form_factor/NormalizedFormFactor.h>
#include <data/Molecule.h>
#include <grid/Grid.h>
#include <dataset/SimpleDataset.h>
#include <plots/All.h>
#include <settings/All.h>

#include "hist/hist_test_helper.h"
#include "grid/grid_debug.h"

using namespace ausaxs;
using namespace ausaxs::data;

TEST_CASE("CompositeDistanceHistogramFFGrid::volumes", "[manual]") {
    settings::general::verbose = false;
    data::Molecule protein("tests/files/2epe.pdb");

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

auto calc_scat = [] (double k) {
    const auto& q_axis = constants::axes::q_vals;
    auto ff_C = form_factor::lookup::atomic::raw::get(form_factor::form_factor_t::C);

    auto V = std::pow(settings::grid::exv::width, 3);
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
    auto ff_C = form_factor::lookup::atomic::raw::get(form_factor::form_factor_t::C);
    auto ff_O = form_factor::lookup::atomic::raw::get(static_cast<form_factor::form_factor_t>(form_factor::water_bin));
    form_factor::FormFactor ffx = form_factor::ExvFormFactor(std::pow(settings::grid::exv::width, 3));
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
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = true;
    settings::hist::weighted_bins = true;
    settings::grid::exv::width = 1;
    settings::grid::min_exv_radius = 0;

    std::vector<AtomFF> atoms = SimpleCube::get_atoms();
    atoms.emplace_back(AtomFF({0, 0, 0}, form_factor::form_factor_t::C));

    std::vector<Water> waters = SimpleCube::get_waters();
    waters.emplace_back(Water({0, 0, 0}));

    Molecule protein({Body{atoms, waters}});
    GridDebug::generate_debug_grid(protein);

    SECTION("Grid") {
        auto h = DebugHistogramManagerMTFFGrid<false>(&protein).calculate_all();
        auto h_cast = static_cast<hist::CompositeDistanceHistogramFFGrid*>(h.get());
        REQUIRE(SimpleCube::check_exact(h_cast->get_d_axis()));
        REQUIRE(SimpleCube::check_exact(h_cast->get_d_axis_ax()));
        REQUIRE(SimpleCube::check_exact(h_cast->get_d_axis_xx()));

        auto Iq_exp = calc_scat_water();
        REQUIRE(compare_hist(Iq_exp, h->debye_transform()));
    }

    SECTION("GridSurface") {
        auto h = DebugHistogramManagerMTFFGridSurface<false>(&protein).calculate_all();
        auto h_cast = static_cast<hist::CompositeDistanceHistogramFFGridSurface*>(h.get());
        REQUIRE(SimpleCube::check_exact(h_cast->get_d_axis()));
        REQUIRE(SimpleCube::check_exact(h_cast->get_d_axis_ax()));
        REQUIRE(SimpleCube::check_exact(h_cast->get_d_axis_xx()));

        auto Iq_exp = calc_scat_water();
        REQUIRE(compare_hist(Iq_exp, h->debye_transform()));
    }

    SECTION("GridScalableExv") {
        auto h = DebugHistogramManagerMTFFGridScalableExv<false>(&protein).calculate_all();
        auto h_cast = static_cast<hist::CompositeDistanceHistogramFFGridScalableExv*>(h.get());
        REQUIRE(SimpleCube::check_exact(h_cast->get_d_axis()));
        REQUIRE(SimpleCube::check_exact(h_cast->get_d_axis_ax()));
        REQUIRE(SimpleCube::check_exact(h_cast->get_d_axis_xx()));

        auto Iq_exp = calc_scat_water();
        REQUIRE(compare_hist(Iq_exp, h->debye_transform()));
    }
}

// Check that solvent scaling is consistent
TEST_CASE("HistogramManagerMTFFGrid: consistent solvent density fitting", "[files]") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = true;
    settings::hist::weighted_bins = true;
    settings::grid::min_exv_radius = 0;
    settings::grid::exv::width = 1;

    data::Molecule protein("tests/files/2epe.pdb");
    auto hg = hist::HistogramManagerMTFFGrid<false>(&protein).calculate_all();
    auto hg_cast = static_cast<hist::CompositeDistanceHistogramFFGrid*>(hg.get());
    
    auto hs = hist::HistogramManagerMTFFGridSurface<false>(&protein).calculate_all();
    auto hs_cast = static_cast<hist::CompositeDistanceHistogramFFGridSurface*>(hs.get());

    auto hse = hist::HistogramManagerMTFFGridScalableExv<false>(&protein).calculate_all();
    auto hse_cast = static_cast<hist::CompositeDistanceHistogramFFGridScalableExv*>(hse.get());

    std::vector<double> rho = {0.2, 0.3, 0.4, 0.5};
    for (auto r : rho) {
        hg_cast->apply_solvent_density_scaling_factor(r);
        hs_cast->apply_solvent_density_scaling_factor(r);
        hse_cast->apply_solvent_density_scaling_factor(r);

        auto Ig  = hg_cast->debye_transform();
        auto Is  = hs_cast->debye_transform();
        auto Ise = hse_cast->debye_transform();

        REQUIRE(compare_hist(Ig, Is));
        REQUIRE(compare_hist(Ig, Ise));
    }
}

// Check that the surface scaling works as expected
TEST_CASE("HistogramManagerMTFFGridSurface: surface_scaling") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = true;
    settings::hist::weighted_bins = true;
    settings::grid::cell_width = 1;
    settings::grid::exv::width = 1;
    settings::grid::min_exv_radius = 0;

    std::vector<AtomFF> atoms = SimpleCube::get_atoms();
    atoms.emplace_back(AtomFF({0, 0, 0}, form_factor::form_factor_t::C));

    Molecule protein({Body{atoms}});
    GridDebug::generate_debug_grid(protein);
    auto h = DebugHistogramManagerMTFFGridSurface<false>(&protein).calculate_all();
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

// Check that the excluded volume scaling works as expected
TEST_CASE("HistogramManagerMTFFGridScalableExv: exv scaling") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = true;
    settings::hist::weighted_bins = true;
    settings::grid::cell_width = 1;
    settings::grid::exv::width = 1;
    settings::grid::min_exv_radius = 0;

    SECTION("simple") {
        std::vector<AtomFF> atoms = {AtomFF({0, 0, 0}, form_factor::form_factor_t::C)};

        Molecule protein({Body{atoms}});
        GridDebug::generate_debug_grid(protein); // overrides exv generation to a known configuration
        auto h = DebugHistogramManagerMTFFGridScalableExv<false>(&protein).calculate_all();
        auto h_cast = static_cast<hist::CompositeDistanceHistogramFFGridScalableExv*>(h.get());

        auto calc = [] (double k) {
            const auto& q_axis = constants::axes::q_vals;
            auto ff_C = form_factor::lookup::atomic::raw::get(form_factor::form_factor_t::C);
            auto ffx = form_factor::ExvFormFactor(std::pow(settings::grid::exv::width*k, 3));
            auto d = SimpleCube::d_exact;
            std::for_each(d.begin(), d.end(), [k] (double& v) {v *= k;});

            std::vector<double> Iq_exp(q_axis.size(), 0);
            for (unsigned int q = 0; q < q_axis.size(); ++q) {
                // calculation ax: 1 x 9 points
                //          1 line  of length 0
                //          8 lines of length sqrt(3) = 1.73
                //
                // calculation xx: 9 x 9 points
                //          1 line  of length 0
                //          3 lines of length 2
                //          3 lines of length sqrt(2*2^2) = sqrt(8) = 2.82
                //          1 line  of length sqrt(3*2^2) = sqrt(12) = 3.46

                double aasum = 
                    1;
                double axsum = 
                    1 +
                    8*std::sin(q_axis[q]*d[1])/(q_axis[q]*d[1]);
                double xxsum = 
                    1 +
                    8 + 
                    16*std::sin(q_axis[q]*d[1])/(q_axis[q]*d[1]) + 
                    24*std::sin(q_axis[q]*d[2])/(q_axis[q]*d[2]) + 
                    24*std::sin(q_axis[q]*d[3])/(q_axis[q]*d[3]) + 
                    8*std::sin(q_axis[q]*d[4])/(q_axis[q]*d[4]);

                Iq_exp[q] += aasum*std::pow(ff_C.evaluate(q_axis[q]), 2);               // + aa
                Iq_exp[q] -= 2*axsum*ff_C.evaluate(q_axis[q])*ffx.evaluate(q_axis[q]);  // -2ax
                Iq_exp[q] += xxsum*std::pow(ffx.evaluate(q_axis[q]), 2);                // + xx
            }
            return Iq_exp;
        };

        // x1
        auto Iq_exp = calc(1);
        h_cast->apply_excluded_volume_scaling_factor(1);
        REQUIRE(compare_hist(Iq_exp, h->debye_transform()));

        // x2
        Iq_exp = calc(2);
        h_cast->apply_excluded_volume_scaling_factor(2);
        REQUIRE(compare_hist(Iq_exp, h->debye_transform()));

        // x3
        Iq_exp = calc(3);
        h_cast->apply_excluded_volume_scaling_factor(3);
        REQUIRE(compare_hist(Iq_exp, h->debye_transform()));

        // x0.5
        Iq_exp = calc(0.5);
        h_cast->apply_excluded_volume_scaling_factor(0.5);
        REQUIRE(compare_hist(Iq_exp, h->debye_transform()));
    }

    SECTION("cube") {
        std::vector<AtomFF> atoms = SimpleCube::get_atoms();
        atoms.emplace_back(AtomFF({0, 0, 0}, form_factor::form_factor_t::C));

        Molecule protein({Body{atoms}});
        GridDebug::generate_debug_grid(protein); // overrides exv generation to a known configuration
        auto h = hist::HistogramManagerMTFFGridScalableExv<false>(&protein).calculate_all();
        auto h_cast = static_cast<hist::CompositeDistanceHistogramFFGridScalableExv*>(h.get());

        auto calc = [] (double k) {
            const auto& q_axis = constants::axes::q_vals;
            auto ff_C = form_factor::lookup::atomic::raw::get(form_factor::form_factor_t::C);
            form_factor::FormFactor ffx = form_factor::ExvFormFactor(std::pow(settings::grid::exv::width*k, 3));
            auto d = SimpleCube::d_exact;

            std::vector<double> Iq_exp(q_axis.size(), 0);
            for (unsigned int q = 0; q < q_axis.size(); ++q) {
                double aasum = 
                    9 + 
                    16*std::sin(q_axis[q]*d[1])/(q_axis[q]*d[1]) + 
                    24*std::sin(q_axis[q]*d[2])/(q_axis[q]*d[2]) + 
                    24*std::sin(q_axis[q]*d[3])/(q_axis[q]*d[3]) + 
                    8 *std::sin(q_axis[q]*d[4])/(q_axis[q]*d[4]);

                double dc = std::sqrt(3);
                double dck = k * std::sqrt(3);
                double d1 = std::sqrt(3) * std::abs(1 - k);
                double d2 = std::sqrt(3 - 2*k + 3*k*k);
                double d3 = std::sqrt(3 + 2*k + 3*k*k);
                double d4 = std::sqrt(3) * (1 + k);
                double axsum = 1 +
                        8 * std::sin(q_axis[q] * dc) / (q_axis[q] * dc) +
                        8 * std::sin(q_axis[q] * dck) / (q_axis[q] * dck) +
                        (k == 1 ? 8 : 8 * std::sin(q_axis[q] * d1) / (q_axis[q] * d1)) +
                        24 * std::sin(q_axis[q] * d2) / (q_axis[q] * d2) +
                        24 * std::sin(q_axis[q] * d3) / (q_axis[q] * d3) +
                        8 * std::sin(q_axis[q] * d4) / (q_axis[q] * d4);

                double xxsum = 
                    9 +
                    16*std::sin(q_axis[q]*d[1]*k)/(q_axis[q]*d[1]*k) + 
                    24*std::sin(q_axis[q]*d[2]*k)/(q_axis[q]*d[2]*k) + 
                    24*std::sin(q_axis[q]*d[3]*k)/(q_axis[q]*d[3]*k) + 
                    8 *std::sin(q_axis[q]*d[4]*k)/(q_axis[q]*d[4]*k);
                
                Iq_exp[q] += aasum*std::pow(ff_C.evaluate(q_axis[q]), 2);               // + aa
                Iq_exp[q] -= 2*axsum*ff_C.evaluate(q_axis[q])*ffx.evaluate(q_axis[q]);  // -2ax
                Iq_exp[q] += xxsum*std::pow(ffx.evaluate(q_axis[q]), 2);                // + xx
            }
            return Iq_exp;
        };

        // x1
        auto Iq_exp = calc(1);
        h_cast->apply_excluded_volume_scaling_factor(1);
        REQUIRE(compare_hist(Iq_exp, h->debye_transform()));

        // x2
        Iq_exp = calc(2);
        h_cast->apply_excluded_volume_scaling_factor(2);
        REQUIRE(compare_hist(Iq_exp, h->debye_transform()));

        // x3
        Iq_exp = calc(3);
        h_cast->apply_excluded_volume_scaling_factor(3);
        REQUIRE(compare_hist(Iq_exp, h->debye_transform()));

        // x0.5
        Iq_exp = calc(0.5);
        h_cast->apply_excluded_volume_scaling_factor(0.5);
        REQUIRE(compare_hist(Iq_exp, h->debye_transform()));
    }
}