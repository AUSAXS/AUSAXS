#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <em/detail/header/MRCHeader.h>
#include <em/manager/ProteinManagerFactory.h>
#include <em/manager/ProteinManager.h>
#include <em/manager/SimpleProteinManager.h>
#include <em/manager/SmartProteinManager.h>
#include <em/ImageStack.h>
#include <em/Image.h>
#include <data/state/StateManager.h>
#include <data/state/Signaller.h>
#include <data/state/BoundSignaller.h>
#include <data/state/UnboundSignaller.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <fitter/Fit.h>
#include <hist/distance_calculator/HistogramManager.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <hist/HistFwd.h>
#include <settings/All.h>

#include <memory>
#include <iostream>

#include <hist/hist_test_helper.h>

using namespace data;

TEST_CASE("partial_histogram_manager_works") {
    settings::molecule::use_effective_charge = false;

    std::vector<Body> bodies(5);
    Molecule protein(bodies);
    auto phm = protein.get_histogram_manager();
    auto manager = phm->get_state_manager();

    manager->reset_to_false();
    phm->get_probe(0)->external_change();
    phm->get_probe(2)->external_change();
    CHECK(manager->get_externally_modified_bodies() == std::vector{true, false, true, false, false});
}

TEST_CASE("protein_manager") {
    settings::molecule::use_effective_charge = false;

    std::vector<Body> bodies(5);
    Molecule protein(bodies);

    auto phm = protein.get_histogram_manager();
    auto manager = phm->get_state_manager();

    manager->reset_to_false();
    CHECK(manager->get_externally_modified_bodies() == std::vector{false, false, false, false, false});
    std::shared_ptr<signaller::Signaller> probe0 = phm->get_probe(0);
    std::shared_ptr<signaller::Signaller> probe2 = phm->get_probe(2);
    probe0->external_change();
    probe2->external_change();
    CHECK(manager->get_externally_modified_bodies() == std::vector{true, false, true, false, false});

    manager->reset_to_false();
    CHECK(manager->get_externally_modified_bodies() == std::vector{false, false, false, false, false});
    Body body;
    body.register_probe(probe0);
    body.changed_external_state();
    CHECK(manager->get_externally_modified_bodies() == std::vector{true, false, false, false, false});

    manager->reset_to_false();
    CHECK(manager->get_externally_modified_bodies() == std::vector{false, false, false, false, false});
    protein.bind_body_signallers();
    protein.get_body(0).changed_external_state();
    protein.get_body(2).changed_external_state();
    CHECK(manager->get_externally_modified_bodies() == std::vector{true, false, true, false, false});

    manager->reset_to_false();
    CHECK(manager->get_externally_modified_bodies() == std::vector{false, false, false, false, false});
    protein.get_body(0) = Body();
    protein.get_body(4) = Body();
    protein.get_body(0).changed_external_state();
    protein.get_body(4).changed_external_state();
    CHECK(manager->get_externally_modified_bodies() == std::vector{true, false, false, false, true});
}

TEST_CASE("em_partial_histogram_manager") {
    settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::HistogramManagerMT; // don't use the phm since it eats too much memory
    settings::molecule::use_effective_charge = false;

    auto compare = [] (std::shared_ptr<em::managers::ProteinManager> manager1, std::shared_ptr<em::managers::ProteinManager> manager2, double cutoff) {
        auto h1 = manager1->get_histogram(cutoff);
        auto h2 = manager2->get_histogram(cutoff);

        auto p1 = h1->get_total_counts();
        auto p2 = h2->get_total_counts();
        if (p1.size() != p2.size()) {
            unsigned int lim = std::min(p1.size(), p2.size());
            for (unsigned int i = 0; i < lim; i++) {
                std::cout << "h1: " << p1[i] << ", h2: " << p2[i] << std::endl;
            }
            std::cout << "Unequal sizes. " << std::endl;
            return false;
        }
        for (unsigned int i = 0; i < p1.size(); i++) {
            if (p1[i] != p2[i]) {
                std::cout << "Failed on index " << i << ". Values: " << p1[i] << ", " << p2[i] << std::endl;
                return false;
            }
        }
        return true;
    };

    SECTION("basic functionality works") {
        settings::em::fixed_weights = GENERATE(true, false);

        std::shared_ptr<em::detail::header::MRCHeader> header;
        {
            em::detail::header::MRCData data;
            data.cella_x = 1, data.cella_y = 1, data.cella_z = 1, data.nz = 1;
            header = std::make_shared<em::detail::header::MRCHeader>(std::move(data));
        }

        Matrix data = Matrix<float>{{1, 2, 3, 4, 5, 6}, {0.5, 1.5, 2.5, 3.5, 4.5, 5.5}};
        em::Image image(data, header.get(), 0);
        em::ImageStack images({image});

        auto manager = em::factory::create_manager(&images);
        manager->set_charge_levels({2, 4, 6, 8});
        auto protein = manager->get_protein(0);

        REQUIRE(protein->size_body() == 5);
        CHECK(protein->get_body(0).size_atom() == 3);
        CHECK(protein->get_body(1).size_atom() == 4);
        CHECK(protein->get_body(2).size_atom() == 4);
        CHECK(protein->get_body(3).size_atom() == 1);
        CHECK(protein->get_body(4).size_atom() == 0);

        protein = manager->get_protein(3);
        REQUIRE(protein->size_body() == 5);
        CHECK(protein->get_body(0).size_atom() == 0);
        CHECK(protein->get_body(1).size_atom() == 2);
        CHECK(protein->get_body(2).size_atom() == 4);
        CHECK(protein->get_body(3).size_atom() == 1);
        CHECK(protein->get_body(4).size_atom() == 0);
    }

    SECTION("comparison with standard approach") {
        SECTION("simple") {
            std::shared_ptr<em::detail::header::MRCHeader> header;
            {
                em::detail::header::MRCData data;
                data.cella_x = 1, data.cella_y = 1, data.cella_z = 1, data.nz = 1;
                header = std::make_shared<em::detail::header::MRCHeader>(std::move(data));
            }

            Matrix data = Matrix<float>{{1, 2, 3, 4, 5, 6}, {0.5, 1.5, 2.5, 3.5, 4.5, 5.5}};
            em::Image image(data, header.get(), 0);
            em::ImageStack images({image});

            std::shared_ptr<em::managers::ProteinManager> manager1 = std::make_shared<em::managers::SimpleProteinManager>(&images);
            std::shared_ptr<em::managers::ProteinManager> manager2 = std::make_shared<em::managers::SmartProteinManager>(&images);
            manager2->set_charge_levels({2, 4, 6, 8});

            // try an arbitrary cutoff level
            REQUIRE(compare(manager1, manager2, 3));

            // try a lower cutoff level
            REQUIRE(compare(manager1, manager2, 1));

            // try a higher cutoff level
            REQUIRE(compare(manager1, manager2, 4));

            // some more tests
            REQUIRE(compare(manager1, manager2, 5));
            REQUIRE(compare(manager1, manager2, 2));
            REQUIRE(compare(manager1, manager2, 3.6));
            REQUIRE(compare(manager1, manager2, 1));
        }

        SECTION("real example") {
            em::ImageStack images("tests/files/A2M_2020_Q4.ccp4");
            std::shared_ptr<em::managers::ProteinManager> manager1 = std::make_shared<em::managers::SimpleProteinManager>(&images);
            std::shared_ptr<em::managers::ProteinManager> manager2 = std::make_shared<em::managers::SmartProteinManager>(&images);

            REQUIRE(compare(manager1, manager2, 4));
            REQUIRE(compare(manager1, manager2, 3));
            REQUIRE(compare(manager1, manager2, 2));
            REQUIRE(compare(manager1, manager2, 5));
            REQUIRE(compare(manager1, manager2, 2));
            REQUIRE(compare(manager1, manager2, 6));
        }
    }
}

#include <data/Molecule.h>
TEST_CASE("SmartProteinManager: consistent_profiles") {
    settings::general::threads = 6;
    settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::HistogramManagerMT;
    em::ImageStack images("tests/files/A2M_2020_Q4.ccp4");

    // SECTION("test protein generator") {
    //     // alpha as the outer loop to ensure the protein is generated anew every time
    //     for (int alpha = 5; alpha < 24; ++alpha) {
    //         images.set_protein_manager(std::make_unique<em::managers::SimpleProteinManager>(&images));
    //         hist::ScatteringProfile hist = images.get_histogram(alpha)->debye_transform();
    //         for (unsigned int charge_levels = 10; charge_levels < 100; charge_levels += 10) {
    //             settings::em::charge_levels = charge_levels;
    //             images.set_protein_manager(std::make_unique<em::managers::SmartProteinManager>(&images));
    //             REQUIRE(images.get_protein_manager()->get_charge_levels().size() == charge_levels+1);
    //             REQUIRE(compare_hist(hist, images.get_histogram(alpha)->debye_transform()));
    //         }
    //     }
    // }

    // SECTION("test protein updater") {
    //     // alpha as the inner loop to check the protein update functionality 
    //     int alpha_min = 5, alpha_max = 24;
    //     std::unordered_map<double, hist::ScatteringProfile> hists;
    //     images.set_protein_manager(std::make_unique<em::managers::SimpleProteinManager>(&images));
    //     for (int alpha = alpha_min; alpha < alpha_max; ++alpha) {
    //         hists[alpha] = images.get_histogram(alpha)->debye_transform();
    //     }

    //     for (unsigned int charge_levels = 10; charge_levels < 100; charge_levels += 10) {
    //         for (int alpha = alpha_min; alpha < alpha_max; ++alpha) {
    //             settings::em::charge_levels = charge_levels;
    //             images.set_protein_manager(std::make_unique<em::managers::SmartProteinManager>(&images));
    //             REQUIRE(images.get_protein_manager()->get_charge_levels().size() == charge_levels+1);
    //             REQUIRE(compare_hist(hists.at(alpha), images.get_histogram(alpha)->debye_transform()));
    //         }
    //     }
    // }

    SECTION("chi2") {
        settings::em::sample_frequency = 2;
        auto res = images.fit("tests/files/A2M_native.dat");
        for (unsigned int charge_levels = 10; charge_levels < 100; charge_levels+= 10) {
            settings::em::charge_levels = charge_levels;
            images.set_protein_manager(std::make_unique<em::managers::SmartProteinManager>(&images));
            REQUIRE(images.get_protein_manager()->get_charge_levels().size() == charge_levels+1);
            REQUIRE_THAT(images.fit("tests/files/A2M_native.dat")->fval, Catch::Matchers::WithinRel(res->fval, 1e-3));
        }
    }
}