#include <catch2/catch_test_macros.hpp>

#include <data/Protein.h>
#include <data/Atom.h>
#include <data/Water.h>
#include <data/Body.h>
#include <data/state/StateManager.h>
#include <data/state/Signaller.h>
#include <hist/HistogramManager.h>
#include <hist/HistogramManagerMT.h>
#include <hist/PartialHistogramManager.h>
#include <hist/PartialHistogramManagerMT.h>
#include <hist/CompositeDistanceHistogramFF.h>
#include <io/ExistingFile.h>
#include <settings/ProteinSettings.h>

struct analytical_histogram {
    static void set_unity_charge(Protein& protein) {
        // set the weights to 1 so we can analytically determine the result
        // waters
        for (auto& atom : protein.get_waters()) {
            atom.set_effective_charge(1);
        }
        // atoms
        for (auto& body : protein.get_bodies()) {
            for (auto& atom : body.get_atoms()) {
                atom.set_effective_charge(1);
            }
        }
    }

    // the following just describes the eight corners of a cube centered at origo
    // {Vector3<double>(-1, -1, -1), Vector3<double>(-1, 1, -1)};
    // {Vector3<double>( 1, -1, -1), Vector3<double>( 1, 1, -1)};
    // {Vector3<double>(-1, -1,  1), Vector3<double>(-1, 1,  1)};
    // {Vector3<double>( 1, -1,  1), Vector3<double>( 1, 1,  1)};

    // calculation: 8 identical points. 
    //      each point has:
    //          1 line  of length 0
    //          3 lines of length 2
    //          3 lines of length sqrt(2*2^2) = sqrt(8) = 2.82
    //          1 line  of length sqrt(3*2^2) = sqrt(12) = 3.46
    std::vector<double> p_exp = {8, 0, 2*8*3, 8, 0, 0, 0, 0, 0, 0};
};

/**
 * @brief Compare two histograms. 
 *        Only indices [0, p1.size()] are checked.
 */
bool compare_hist(Vector<double> p1, Vector<double> p2) {
    for (unsigned int i = 0; i < p1.size(); i++) {
        if (!utility::approx(p1[i], p2[i])) {
            std::cout << "Failed on index " << i << ". Values: " << p1[i] << ", " << p2[i] << std::endl;
            return false;
        }
    }
    return true;
}

TEST_CASE_METHOD(analytical_histogram, "HistogramManager::calculate_all") {
    settings::protein::use_effective_charge = false;

    SECTION("analytical") {
        SECTION("atoms only") {
            // the following just describes the eight corners of a cube centered at origo
            std::vector<Atom> b1 = {Atom(Vector3<double>(-1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, "C", "C", 1)};
            std::vector<Atom> b2 = {Atom(Vector3<double>( 1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>( 1, 1, -1), 1, "C", "C", 1)};
            std::vector<Atom> b3 = {Atom(Vector3<double>(-1, -1,  1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1,  1), 1, "C", "C", 1)};
            std::vector<Atom> b4 = {Atom(Vector3<double>( 1, -1,  1), 1, "C", "C", 1), Atom(Vector3<double>( 1, 1,  1), 1, "C", "C", 1)};
            std::vector<Body> a = {Body(b1), Body(b2), Body(b3), Body(b4)};
            Protein protein(a);
            set_unity_charge(protein);

            { // hm
                auto hm = hist::HistogramManager(&protein).calculate_all();
                REQUIRE(compare_hist(p_exp, hm->p));
            }
            { // hm_mt
                auto hm_mt = hist::HistogramManagerMT(&protein).calculate_all();
                REQUIRE(compare_hist(p_exp, hm_mt->p));
            }
            { // phm
                auto phm = protein.get_histogram();
                REQUIRE(compare_hist(p_exp, phm->p));
            }
            { // phm_mt
                auto phm_mt = hist::PartialHistogramManagerMT(&protein).calculate_all();
                REQUIRE(compare_hist(p_exp, phm_mt->p));
            }
        }

        SECTION("waters only") {
            // the following just describes the eight corners of a cube centered at origo, with an additional atom at the very middle
            std::vector<Atom> a = {};
            std::vector<Water> w = {Water(Vector3<double>(-1, -1, -1), 1, "C", "C", 1), Water(Vector3<double>(-1, 1, -1), 1, "C", "C", 1), 
                                    Water(Vector3<double>( 1, -1, -1), 1, "C", "C", 1), Water(Vector3<double>( 1, 1, -1), 1, "C", "C", 1), 
                                    Water(Vector3<double>(-1, -1,  1), 1, "C", "C", 1), Water(Vector3<double>(-1, 1,  1), 1, "C", "C", 1),
                                    Water(Vector3<double>( 1, -1,  1), 1, "C", "C", 1), Water(Vector3<double>( 1, 1,  1), 1, "C", "C", 1)};
            Protein protein(a, w);
            set_unity_charge(protein);

            { // hm
                auto hm = hist::HistogramManager(&protein).calculate_all();
                REQUIRE(compare_hist(p_exp, hm->p));
            }
            { // hm_mt
                auto hm_mt = hist::HistogramManagerMT(&protein).calculate_all();
                REQUIRE(compare_hist(p_exp, hm_mt->p));
            }
            { // phm
                auto phm = protein.get_histogram();
                REQUIRE(compare_hist(p_exp, phm->p));
            }
            { // phm_mt
                auto phm_mt = hist::PartialHistogramManagerMT(&protein).calculate_all();
                REQUIRE(compare_hist(p_exp, phm_mt->p));
            }
        }

        SECTION("both waters and atoms") {
            // the following just describes the eight corners of a cube centered at origo, with an additional atom at the very middle
            std::vector<Atom> b1 = {Atom( Vector3<double>(-1, -1, -1), 1, "C", "C", 1), Atom( Vector3<double>(-1, 1, -1), 1, "C", "C", 1)};
            std::vector<Atom> b2 = {Atom( Vector3<double>( 1, -1, -1), 1, "C", "C", 1), Atom( Vector3<double>( 1, 1, -1), 1, "C", "C", 1)};
            std::vector<Atom> b3 = {Atom( Vector3<double>(-1, -1,  1), 1, "C", "C", 1), Atom( Vector3<double>(-1, 1,  1), 1, "C", "C", 1)};
            std::vector<Water> w = {Water(Vector3<double>( 1, -1,  1), 1, "C", "C", 1), Water(Vector3<double>( 1, 1,  1), 1, "C", "C", 1)};
            std::vector<Body> a = {Body(b1), Body(b2), Body(b3)};
            Protein protein(a, w);
            set_unity_charge(protein);

            { // hm
                auto hm = hist::HistogramManager(&protein).calculate_all();
                REQUIRE(compare_hist(p_exp, hm->p));
            }
            { // hm_mt
                auto hm_mt = hist::HistogramManagerMT(&protein).calculate_all();
                REQUIRE(compare_hist(p_exp, hm_mt->p));
            }
            { // phm
                auto phm = protein.get_histogram();
                REQUIRE(compare_hist(p_exp, phm->p));
            }
            { // phm_mt
                auto phm_mt = hist::PartialHistogramManagerMT(&protein).calculate_all();
                REQUIRE(compare_hist(p_exp, phm_mt->p));
            }
        }
    }

    SECTION("real data") {
        settings::protein::use_effective_charge = true;
                
        // create the atom, and perform a sanity check on our extracted list
        Protein protein("test/files/2epe.pdb");
        protein.generate_new_hydration();

        auto p_exp = protein.get_histogram()->p;

        { // hm
            auto hm = hist::HistogramManager(&protein).calculate_all();
            REQUIRE(p_exp.size() == hm->p.size());
            REQUIRE(compare_hist(p_exp, hm->p));
        }
        { // hm_mt
            auto hm_mt = hist::HistogramManagerMT(&protein).calculate_all();
            REQUIRE(p_exp.size() == hm_mt->p.size());
            REQUIRE(compare_hist(p_exp, hm_mt->p));
        }
        { // phm
            auto phm = hist::PartialHistogramManager(&protein).calculate_all();
            REQUIRE(p_exp.size() == phm->p.size());
            REQUIRE(compare_hist(p_exp, phm->p));
        }
        { // phm_mt
            auto phm_mt = hist::PartialHistogramManagerMT(&protein).calculate_all();
            REQUIRE(p_exp.size() == phm_mt->p.size());
            REQUIRE(compare_hist(p_exp, phm_mt->p));
        }
    }
}

TEST_CASE("PartialHistogramManager::get_probe") {
    Protein protein("test/files/2epe.pdb");
    auto phm = hist::PartialHistogramManager(&protein);
    auto sm = phm.get_state_manager();

    // check the signalling object is correct
    CHECK(phm.get_probe(0) == sm->get_probe(0)); 

    // check that it links to the state manager
    sm->reset();
    phm.get_probe(0)->external_change();
    CHECK(sm->is_externally_modified(0));
}

TEST_CASE("PartialHistogramManager::signal_modified_hydration_layer") {
    Protein protein("test/files/2epe.pdb");
    auto phm = hist::PartialHistogramManager(&protein);
    auto sm = phm.get_state_manager();
    sm->reset();
    phm.signal_modified_hydration_layer();
    CHECK(sm->get_modified_hydration());
}