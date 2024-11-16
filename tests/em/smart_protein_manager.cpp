#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <hydrate/generation/RadialHydration.h>
#include <em/manager/SimpleProteinManager.h>
#include <em/manager/SmartProteinManager.h>
#include <em/ImageStack.h>
#include <data/Molecule.h>
#include <settings/All.h>

#include <hist/hist_test_helper.h>

using namespace ausaxs;

struct fixture {
    fixture() {
        settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::HistogramManagerMT;
        stack = std::make_unique<em::ImageStack>("tests/files/A2M_2020_Q4.ccp4");
        manager = std::make_unique<em::managers::SmartProteinManager>(stack.get());
    }

    std::unique_ptr<em::ImageStack> stack;
    std::unique_ptr<em::managers::SmartProteinManager> manager;
};

TEST_CASE_METHOD(fixture, "SmartProteinManager::SmartProteinManager", "[files]") {
    SECTION("ImageStackBase*") {
        CHECK(!manager->get_charge_levels().empty());
    }
}

TEST_CASE_METHOD(fixture, "SmartProteinManager::set_charge_levels") {
    std::vector<double> charge_levels = {1, 2, 3};
    manager->set_charge_levels(charge_levels);
    CHECK(manager->get_charge_levels() == std::vector<double>{1, 2, 3, 10000});
}

TEST_CASE_METHOD(fixture, "SmartProteinManager::get_protein", "[files]") {
    // we just check the size of the returned protein
    unsigned int size = manager->get_protein(2)->size_atom();
    CHECK(size != 0);
    CHECK(manager->get_protein(1)->size_atom() > size);
    CHECK(manager->get_protein(3)->size_atom() < size);
    CHECK(manager->get_protein(2)->size_atom() == size);
}

TEST_CASE_METHOD(fixture, "SmartProteinManager::get_histogram", "[files]") {
    CHECK(manager->get_histogram(1)->get_total_counts() == manager->get_protein(1)->get_histogram()->get_total_counts());
}

TEST_CASE("SmartProteinManager::generate_protein", "[files]") {
    settings::general::threads = 6;
    settings::em::sample_frequency = 2;
    settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::PartialHistogramManagerMT;

    // alpha as the outer loop to ensure the protein is generated anew every time
    em::ImageStack images("tests/files/A2M_2020_Q4.ccp4");
    for (int alpha = 5; alpha < 24; ++alpha) {
        images.set_protein_manager(std::make_unique<em::managers::SimpleProteinManager>(&images));
        hist::ScatteringProfile hist = images.get_histogram(alpha)->debye_transform();
        for (unsigned int charge_levels = 10; charge_levels < 100; charge_levels += 10) {
            settings::em::charge_levels = charge_levels;
            images.set_protein_manager(std::make_unique<em::managers::SmartProteinManager>(&images));
            REQUIRE(images.get_protein_manager()->get_charge_levels().size() == charge_levels+1);
            REQUIRE(compare_hist(hist, images.get_histogram(alpha)->debye_transform()));
        }
    }
}

TEST_CASE("SmartProteinManager::update_protein", "[files]") {
    settings::general::threads = 6;
    settings::em::sample_frequency = 2;
    settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::PartialHistogramManagerMT;

    // alpha as the inner loop to check the protein update functionality 
    int alpha_min = 5, alpha_max = 24;
    em::ImageStack images("tests/files/A2M_2020_Q4.ccp4");
    std::unordered_map<double, hist::ScatteringProfile> hists;
    images.set_protein_manager(std::make_unique<em::managers::SimpleProteinManager>(&images));
    for (int alpha = alpha_min; alpha < alpha_max; ++alpha) {
        hists[alpha] = images.get_histogram(alpha)->debye_transform();
    }

    for (unsigned int charge_levels = 10; charge_levels < 50; charge_levels += 10) {
        for (int alpha = alpha_min; alpha < alpha_max; ++alpha) {
            settings::em::charge_levels = charge_levels;
            images.set_protein_manager(std::make_unique<em::managers::SmartProteinManager>(&images));
            REQUIRE(images.get_protein_manager()->get_charge_levels().size() == charge_levels+1);
            REQUIRE(compare_hist(hists.at(alpha), images.get_histogram(alpha)->debye_transform()));
        }
    }
}

TEST_CASE("SmartProteinManager: consistency", "[files]") {
    settings::general::threads = 6;
    settings::general::verbose = false;
    settings::general::supplementary_plots = false;
    settings::em::sample_frequency = 2;
    settings::em::alpha_levels = {6, 8};
    settings::em::save_pdb = false;
    settings::fit::max_iterations = 20;
    settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::PartialHistogramManagerMT;
    hydrate::RadialHydration::set_noise_generator([] () {return Vector3<double>{0, 0, 0};});

    // we need a hydration-sensitive map for this test
    em::ImageStack images("tests/files/emd_24889.map");
    auto res = images.fit("tests/files/SASDJG5.dat");
    for (unsigned int charge_levels = 10; charge_levels < 50; charge_levels+= 10) {
        settings::em::charge_levels = charge_levels;
        images.set_protein_manager(std::make_unique<em::managers::SmartProteinManager>(&images));
        REQUIRE(images.get_protein_manager()->get_charge_levels().size() == charge_levels+1);
        REQUIRE_THAT(images.fit("tests/files/SASDJG5.dat")->fval, Catch::Matchers::WithinRel(res->fval, 1e-3));
    }
}