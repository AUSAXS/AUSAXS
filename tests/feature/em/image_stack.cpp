#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <em/ImageStack.h>
#include <em/detail/header/MRCHeader.h>
#include <em/manager/SmartProteinManager.h>
#include <settings/All.h>
#include <data/Molecule.h>
#include <grid/Grid.h>
#include <hist/histogram_manager/HistogramManagerMT.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <plots/All.h>

#include <map>

using namespace ausaxs;

// TODO: fix this test; it's a solid one
TEST_CASE("ImageStack: test with sphere", "[broken]") {
    settings::general::verbose = false;
    settings::molecule::center = false;

    // generate big sphere
    auto lims = Limit3D(-50, 50, -50, 50, -50, 50);
    grid::Grid grid(lims);
    double radius = 15;
    double radius2 = radius*radius;
    auto axes = grid.get_axes();
    Vector3<double> center = grid.to_xyz(grid.get_center());
    for (unsigned int i = 0; i < axes.x.bins; ++i) {
        for (unsigned int j = 0; j < axes.y.bins; ++j) {
            for (unsigned int k = 0; k < axes.z.bins; ++k) {
                if (grid.to_xyz(i, j, k).distance2(center) < radius2) {
                    grid.grid.index(i, j, k) = grid::detail::VOLUME;
                }
            }
        }
    }
    auto loc = "temp/tests/em/sphere.pdb";
    grid.save(loc);

    data::Molecule protein(loc);
    auto Iq = hist::HistogramManagerMT<true>(&protein).calculate_all()->debye_transform();
    Iq.as_dataset().save("temp/tests/em/sphere_Iq.dat");

    std::unique_ptr<em::detail::header::MRCHeader> header = std::make_unique<em::detail::header::MRCHeader>();
    std::unique_ptr<em::detail::header::MRCData> header_data = std::make_unique<em::detail::header::MRCData>();
    header_data->cella_x = axes.x.span();
    header_data->cella_y = axes.y.span();
    header_data->cella_z = axes.z.span();
    header_data->nx = axes.x.bins;
    header_data->ny = axes.y.bins;
    header_data->nz = axes.z.bins;

    std::vector<em::Image> images(lims.z.span()/settings::grid::cell_width, Matrix<float>(0, 0));
    for (unsigned int k = 0; k < images.size(); ++k) {
        Matrix<float> data(axes.x.bins, axes.y.bins);
        for (unsigned int i = 0; i < axes.x.bins; ++i) {
            for (unsigned int j = 0; j < axes.y.bins; ++j) {
                double dist = std::sqrt(grid.to_xyz(i, j, k).distance2(center));
                data.index(i, j) = radius/dist;
            }
        }
        images[k] = em::Image(data, header.get(), k);
    }

    em::ImageStack stack(images);
    // auto[fit, landscape] = stack.cutoff_scan_fit(100, std::move(Iq));
    auto fit = stack.fit("temp/tests/em/sphere_Iq.dat");
    // plots::PlotLandscape::quick_plot(landscape, "temp/tests/em/sphere_landscape.png");
    REQUIRE(fit->fval/fit->dof < 1.1);
}

TEST_CASE("ImageStack::get_protein") {
    settings::molecule::center = false;
    settings::grid::min_bins = 100;
    
    SECTION("em_weights") {
        SECTION("dynamic") {
            settings::em::fixed_weights = false;
            auto data  = Matrix<float>{{0, 1, 3, 5, 1, 0}, {0, 3, 5, 5, 0, 0}, {0, 0, 1, 3, 3, 0}, {0, 3, 0, 5, 1, 0}, {0, 1, 3, 5, 0, 0}, {0, 1, 0, 3, 1, 5}};
            auto data2 = Matrix<float>{{0, 1, 2, 3, 2, 1}, {0, 3, 2, 1, 3, 0}, {0, 1, 2, 0, 1, 0}, {2, 0, 0, 3, 1, 0}, {0, 1, 2, 1, 1, 0}, {3, 3, 3, 2, 1, 1}};
            em::ImageStack images({data, data2});

            std::unique_ptr header = std::make_unique<em::detail::header::MRCHeader>();
            std::unique_ptr header_data = std::make_unique<em::detail::header::MRCData>();
            header_data->cella_x = 6; header_data->cella_y = 6; header_data->cella_z = 2;
            header_data->nx = 6; header_data->ny = 6; header_data->nz = 2;
            header->set_data(std::move(header_data));
            images.set_header(std::move(header));

            auto protein = images.get_protein(1);
            REQUIRE(protein->size_atom() == 4+3+3+3+3+4 + 5+4+3+3+4+6);
            std::map<float, unsigned int> counts = {{1, 0}, {2, 0}, {3, 0}, {4, 0}, {5, 0}};
            for (const auto& atom : protein->get_atoms()) {
                ++counts[atom.weight()];
            }
            REQUIRE(counts.at(1) == 2+0+1+1+1+2 + 2+1+2+1+3+2);
            REQUIRE(counts.at(2) == 0+0+0+0+0+0 + 2+1+1+1+1+1);
            REQUIRE(counts.at(3) == 1+1+2+1+1+1 + 1+2+0+1+0+3);
            REQUIRE(counts.at(4) == 0);
            REQUIRE(counts.at(5) == 1+2+0+1+1+1 + 0+0+0+0+0+0);
        }

        SECTION("fixed") {
            settings::em::fixed_weights = true;
            auto data  = Matrix<float>{{0, 1, 3, 5, 1, 0}, {0, 3, 5, 5, 0, 0}, {0, 0, 1, 3, 3, 0}, {0, 3, 0, 5, 1, 0}, {0, 1, 3, 5, 0, 0}, {0, 1, 0, 3, 1, 5}};
            auto data2 = Matrix<float>{{0, 1, 2, 3, 2, 1}, {0, 3, 2, 1, 3, 0}, {0, 1, 2, 0, 1, 0}, {2, 0, 0, 3, 1, 0}, {0, 1, 2, 1, 1, 0}, {3, 3, 3, 2, 1, 1}};
            em::ImageStack images({data, data2});

            std::unique_ptr header = std::make_unique<em::detail::header::MRCHeader>();
            std::unique_ptr header_data = std::make_unique<em::detail::header::MRCData>();
            header_data->cella_x = 6; header_data->cella_y = 6; header_data->cella_z = 2;
            header_data->nx = 6; header_data->ny = 6; header_data->nz = 2;
            header->set_data(std::move(header_data));
            images.set_header(std::move(header));
            
            auto protein = images.get_protein(1);
            REQUIRE(protein->get_atoms().size() == 4+3+3+3+3+4 + 5+4+3+3+4+6);
            for (const auto& atom : protein->get_atoms()) {
                REQUIRE(atom.weight() == 1);
            }
            settings::em::fixed_weights = false;
        }
    }
}

// Check that the mass is calculated correctly
TEST_CASE("ImageStack::get_mass") {
    settings::em::sample_frequency = 2;
    settings::em::alpha_levels = GENERATE(Limit{3, 6}, Limit{6, 8}, Limit{8, 14});

    em::ImageStack images("tests/files/A2M_2020_Q4.ccp4");
    std::unordered_map<double, double> vals;
    for (int i = 5; i < 12; ++i) {vals[i] = images.get_mass(images.from_level(i));}
    for (unsigned int charge_levels = 10; charge_levels < 50; charge_levels += 10) {
        settings::em::charge_levels = charge_levels;
        images.set_protein_manager(std::make_unique<em::managers::SmartProteinManager>(&images));

        for (int i = 5; i < 12; ++i) {
            REQUIRE_THAT(images.get_mass(images.from_level(i)), Catch::Matchers::WithinRel(vals.at(i), 1e-3));
        }

        for (int i = 12; i < 5; --i) {
            REQUIRE_THAT(images.get_mass(images.from_level(i)), Catch::Matchers::WithinRel(vals.at(i), 1e-3));
        }
    }
}