#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <em/ImageStack.h>
#include <em/detail/header/MRCHeader.h>
#include <settings/All.h>
#include <data/Atom.h>
#include <data/Protein.h>

#include <map>

TEST_CASE("ImageStack::ImageStack") {
    SECTION("ExistingFile") {}
    SECTION("std::vector<Image>&") {}
    CHECK(false);
}

TEST_CASE("ImageStack::fit") {
    CHECK(false);
}

TEST_CASE("ImageStack::cutoff_scan") {
    CHECK(false);
}

TEST_CASE("ImageStack::cutoff_scan_fit") {
    CHECK(false);
}

TEST_CASE("ImageStack::get_fitted_water_factors") {
    CHECK(false);
}

TEST_CASE("ImageStack::get_fitted_water_factors_dataset") {
    CHECK(false);
}

TEST_CASE("ImageStack::get_protein") {
    SECTION("em_weights") {
        SECTION("dynamic") {
            Matrix data = Matrix<float>{{0, 1, 3, 5, 1, 0}, {0, 3, 5, 5, 0, 0}, {0, 0, 1, 3, 3, 0}, {0, 3, 0, 5, 1, 0}, {0, 1, 3, 5, 0, 0}, {0, 1, 0, 3, 1, 5}};
            Matrix data2 = Matrix<float>{{0, 1, 2, 3, 2, 1}, {0, 3, 2, 1, 3, 0}, {0, 1, 2, 0, 1, 0}, {2, 0, 0, 3, 1, 0}, {0, 1, 2, 1, 1, 0}, {3, 3, 3, 2, 1, 1}};
            em::ImageStack images({data, data2});

            std::unique_ptr header = std::make_unique<em::detail::header::MRCHeader>();
            std::unique_ptr header_data = std::make_unique<em::detail::header::MRCData>();
            header_data->cella_x = 6; header_data->cella_y = 6; header_data->cella_z = 2;
            header->set_data(std::move(header_data));
            images.set_header(std::move(header));

            auto protein = images.get_protein(1);
            REQUIRE(protein->get_atoms().size() == 4+3+3+3+3+4 + 5+4+3+3+4+6);
            std::map<float, unsigned int> counts = {{1, 0}, {2, 0}, {3, 0}, {4, 0}, {5, 0}};
            for (const auto& atom : protein->get_atoms()) {
                counts[atom.get_occupancy()]++;
            }
            REQUIRE(counts.at(1) == 2+0+1+1+1+2 + 2+1+2+1+3+2);
            REQUIRE(counts.at(2) == 0+0+0+0+0+0 + 2+1+1+1+1+1);
            REQUIRE(counts.at(3) == 1+1+2+1+1+1 + 1+2+0+1+0+3);
            REQUIRE(counts.at(4) == 0);
            REQUIRE(counts.at(5) == 1+2+0+1+1+1 + 0+0+0+0+0+0);
        }

        SECTION("fixed") {
            settings::em::fixed_weights = true;
            Matrix data = Matrix<float>{{0, 1, 3, 5, 1, 0}, {0, 3, 5, 5, 0, 0}, {0, 0, 1, 3, 3, 0}, {0, 3, 0, 5, 1, 0}, {0, 1, 3, 5, 0, 0}, {0, 1, 0, 3, 1, 5}};
            Matrix data2 = Matrix<float>{{0, 1, 2, 3, 2, 1}, {0, 3, 2, 1, 3, 0}, {0, 1, 2, 0, 1, 0}, {2, 0, 0, 3, 1, 0}, {0, 1, 2, 1, 1, 0}, {3, 3, 3, 2, 1, 1}};
            em::ImageStack images({data, data2});

            std::unique_ptr header = std::make_unique<em::detail::header::MRCHeader>();
            std::unique_ptr header_data = std::make_unique<em::detail::header::MRCData>();
            header_data->cella_x = 6; header_data->cella_y = 6; header_data->cella_z = 2;
            header->set_data(std::move(header_data));
            images.set_header(std::move(header));
            
            auto protein = images.get_protein(1);
            REQUIRE(protein->get_atoms().size() == 4+3+3+3+3+4 + 5+4+3+3+4+6);
            for (const auto& atom : protein->get_atoms()) {
                REQUIRE(atom.get_occupancy() == 1);
            }
            settings::em::fixed_weights = false;
        }
    }
}