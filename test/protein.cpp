#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

#include <data/Protein.h>
#include <hydrate/Grid.h>
#include <utility/Constants.h>
#include <utility/Utility.h>
#include <fitter/SimpleIntensityFitter.h>
#include <hist/HistogramManagerMT.h>
#include <hist/PartialHistogramManagerMT.h>
#include <plots/all.h>

using std::cout, std::endl, std::vector, std::shared_ptr;

TEST_CASE("simulate_dataset", "[protein],[files]") {
    setting::axes::qmax = 0.4;
    setting::protein::use_effective_charge = false;
    setting::em::sample_frequency = 2;
    Protein protein("test/files/2epe.pdb");
    SimpleDataset data = protein.simulate_dataset();

    SimpleIntensityFitter fitter(data, protein.get_histogram());
    auto res = fitter.fit();
    REQUIRE_THAT(res->fval/res->dof, Catch::Matchers::WithinAbs(1., 0.5));
    plots::PlotIntensityFit plot1(res);
    plot1.save("figures/test/protein/check_chi2_1.png");
}

TEST_CASE("compare_debye", "[protein]") {
    vector<Atom> atoms = {Atom(Vector3<double>(-1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, "C", "C", 1),
                       Atom(Vector3<double>(1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>(1, 1, -1), 1, "C", "C", 1),
                       Atom(Vector3<double>(-1, -1, 1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1, 1), 1, "C", "C", 1),
                       Atom(Vector3<double>(1, -1, 1), 1, "C", "C", 1), Atom(Vector3<double>(1, 1, 1), 1, "C", "C", 1)};
    Protein protein(atoms, {});

    vector<double> I_dumb = protein.calc_debye_scattering_intensity();
    vector<double> I_smart = protein.get_histogram().calc_debye_scattering_intensity().col("I");

    for (int i = 0; i < 8; i++) {
        if (!utility::approx(I_dumb[i], I_smart[i], 1e-1)) {
            cout << "Failed on index " << i << ". Values: " << I_dumb[i] << ", " << I_smart[i] << endl;
            REQUIRE(false);
        }
    }
}

TEST_CASE("compare_debye_real", "[protein],[files],[slow]") {
    Protein protein("test/files/2epe.pdb");
    protein.clear_hydration();

    std::cout << "hydration atoms: " << protein.hydration_atoms.size() << std::endl; 

    vector<double> I_dumb = protein.calc_debye_scattering_intensity();
    vector<double> I_smart = protein.get_histogram().calc_debye_scattering_intensity().col("I");

    for (int i = 0; i < 8; i++) {
        if (!utility::approx(I_dumb[i], I_smart[i], 1e-3, 0.05)) {
            cout << "Failed on index " << i << ". Values: " << I_dumb[i] << ", " << I_smart[i] << endl;
            REQUIRE(false);
        }
    }
}

/**
 * @brief Compare two histograms. 
 *        Only indices [0, p1.size()] are checked.
 */
bool compare_hist(Vector<double> p1, Vector<double> p2) {
    for (size_t i = 0; i < p1.size(); i++) {
        if (!utility::approx(p1[i], p2[i])) {
            cout << "Failed on index " << i << ". Values: " << p1[i] << ", " << p2[i] << endl;
            return false;
        }
    }
    return true;
}

TEST_CASE("histogram", "[protein]") {
    setting::axes::scattering_intensity_plot_binned_width = 1;

    SECTION("multiple bodies, simple") {
        setting::protein::use_effective_charge = true;
        // make the protein
        vector<Atom> b1 = {Atom(Vector3<double>(-1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, "C", "C", 1)};
        vector<Atom> b2 = {Atom(Vector3<double>(1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>(1, 1, -1), 1, "C", "C", 1)};
        vector<Atom> b3 = {Atom(Vector3<double>(-1, -1, 1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1, 1), 1, "C", "C", 1)};
        vector<Atom> b4 = {Atom(Vector3<double>(1, -1, 1), 1, "C", "C", 1), Atom(Vector3<double>(1, 1, 1), 1, "C", "C", 1)};
        vector<vector<Atom>> ap = {b1, b2, b3, b4};
        Protein many(ap, {});

        // make the body
        vector<Atom> ab = {Atom(Vector3<double>(-1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, "C", "C", 1),
                            Atom(Vector3<double>(1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>(1, 1, -1), 1, "C", "C", 1),
                            Atom(Vector3<double>(-1, -1, 1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1, 1), 1, "C", "C", 1),
                            Atom(Vector3<double>(1, -1, 1), 1, "C", "C", 1), Atom(Vector3<double>(1, 1, 1), 1, "C", "C", 1)};
        Protein one(ab, {});

        // create some water molecules
        vector<Water> ws(10);
        for (size_t i = 0; i < ws.size(); i++) {
            ws[i] = Water::create_new_water(Vector3<double>(i, i, i));
        }

        many.waters() = ws;
        one.waters() = ws;

        // we now have a protein consisting of three bodies with the exact same contents as a single body.
        // the idea is now to compare the ScatteringHistogram output from their distance calculations, since it
        // is far easier to do for the single body. 
        hist::ScatteringHistogram d_m = many.get_histogram();
        hist::ScatteringHistogram d_o = one.get_histogram();

        // direct access to the histogram data (only p is defined)
        const vector<double>& p_m = d_m.p;
        const vector<double>& p_o = d_o.p;

        // compare each entry
        for (size_t i = 0; i < p_o.size(); i++) {
            if (!utility::approx(p_o[i], p_m[i])) {
                cout << "Failed on index " << i << ". Values: " << p_m[i] << ", " << p_o[i] << endl;
                REQUIRE(false);
            }
        }
        REQUIRE(true);
    }

    SECTION("multiple bodies, real input") {
        setting::protein::use_effective_charge = true;
        Body body("data/lysozyme/2epe.pdb");
        body.center();
        
        // We iterate through the protein data from the body, and split it into multiple pieces of size 100.  
        vector<vector<Atom>> patoms; // vector containing the pieces we split it into
        vector<Atom> p_current(100); // vector containing the current piece
        unsigned int index = 0;      // current index in p_current
        for (unsigned int i = 0; i < body.atoms().size(); i++) {
            p_current[index] = body.atoms(i);
            index++;
            if (index == 100) { // if index is 100, reset to 0
                patoms.push_back(p_current);
                index = 0;
            }
        }

        // add the final few atoms to our list
        if (index != 0) {
            p_current.resize(index);
            patoms.push_back(p_current);
        }

        // create the atom, and perform a sanity check on our extracted list
        Protein protein(patoms, {});
        vector<Atom> protein_atoms = protein.atoms();
        vector<Atom> body_atoms = body.atoms();

        // sizes must be equal. this also serves as a separate consistency check on the body generation. 
        if (protein_atoms.size() != body_atoms.size()) {
            cout << "Sizes " << protein_atoms.size() << " and " << body_atoms.size() << " should be equal. " << endl;
            REQUIRE(false);
        }

        // stronger consistency check - we check that all atoms are equal, and appear in the exact same order
        for (unsigned int i = 0; i < protein_atoms.size(); i++) {
            if (protein_atoms[i] != body_atoms[i]) {
                cout << "Comparison failed on index " << i << endl;
                cout << protein_atoms[i].as_pdb() << endl;
                cout << body_atoms[i].as_pdb() << endl;
                REQUIRE(false);
            }
        }

        // generate a hydration layer for the protein, and copy it over to the body
        protein.generate_new_hydration();

        // generate the distance histograms
        hist::ScatteringHistogram d_p = protein.get_histogram();
        hist::ScatteringHistogram d_b = hist::HistogramManager(&protein).calculate_all();

        // direct access to the histogram data (only p is defined)
        const vector<double>& p = d_p.p;
        const vector<double>& b_tot = d_b.p;

        // compare each entry
        for (unsigned int i = 0; i < b_tot.size(); i++) {
            if (!utility::approx(p[i], b_tot[i])) {
                cout << "Failed on index " << i << ". Values: " << p[i] << ", " << b_tot[i] << endl;
                REQUIRE(false);
            }
        }
        REQUIRE(true);
    }

    SECTION("equivalent to old approach") {
        setting::protein::use_effective_charge = false;
        vector<Atom> atoms = {Atom(Vector3<double>(-1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, "C", "C", 1),
                              Atom(Vector3<double>(1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>(1, 1, -1), 1, "C", "C", 1),
                              Atom(Vector3<double>(-1, -1, 1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1, 1), 1, "C", "C", 1),
                              Atom(Vector3<double>(1, -1, 1), 1, "C", "C", 1), Atom(Vector3<double>(1, 1, 1), 1, "C", "C", 1)};

        // new auto-scaling approach
        Protein protein1(atoms);
        Grid grid1(atoms);
        protein1.set_grid(grid1);

        // old approach
        Protein protein2(atoms);
        Axis3D axes(setting::grid::axes, setting::grid::width);
        Grid grid2(axes, setting::grid::width); 
        grid2.add(atoms);
        protein2.set_grid(grid2);

        // generate the distance histograms
        hist::ScatteringHistogram h1 = protein1.get_histogram();
        hist::ScatteringHistogram h2 = protein2.get_histogram();

        // direct access to the histogram data (only p is defined)
        const vector<double>& p1 = h1.p;
        const vector<double>& p2 = h2.p;

        // compare each entry
        for (size_t i = 0; i < p1.size(); i++) {
            if (!utility::approx(p1[i], p2[i])) {
                cout << "Failed on index " << i << ". Values: " << p1[i] << ", " << p2[i] << endl;
                REQUIRE(false);
            }
        }
        REQUIRE(true);
    }
}

TEST_CASE("distance_histograms", "[protein]") {
    setting::protein::use_effective_charge = false;

    SECTION("analytical") {
        SECTION("atoms only") {
            // the following just describes the eight corners of a cube centered at origo, with an additional atom at the very middle
            vector<Atom> b1 = {Atom(Vector3<double>(-1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, "C", "C", 1)};
            vector<Atom> b2 = {Atom(Vector3<double>(1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>(1, 1, -1), 1, "C", "C", 1)};
            vector<Atom> b3 = {Atom(Vector3<double>(-1, -1, 1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1, 1), 1, "C", "C", 1)};
            vector<Water> w = {Water(Vector3<double>(1, -1, 1), 1, "C", "C", 1),   Water(Vector3<double>(1, 1, 1), 1, "C", "C", 1)};
            vector<vector<Atom>> a = {b1, b2, b3};
            Protein protein(a, w);

            // set the weights to 1 so we can analytically determine the result
            // waters
            for (auto& atom : protein.hydration_atoms) {
                atom.set_effective_charge(1);
            }
            // atoms
            for (auto& body : protein.bodies) {
                for (auto& atom : body.atoms()) {
                    atom.set_effective_charge(1);
                }
            }
            protein.updated_charge = true;

            // calculation: 8 identical points. 
            //      each point has:
            //          1 line  of length 0
            //          3 lines of length 2
            //          3 lines of length sqrt(2*2^2) = sqrt(8) = 2.82
            //          1 line  of length sqrt(3*2^2) = sqrt(12) = 3.46
            const vector<double> p_exp = {8, 0, 2*8*3, 8, 0, 0, 0, 0, 0, 0};

            { // hm
                hist::ScatteringHistogram hm = hist::HistogramManager(&protein).calculate_all();
                REQUIRE(compare_hist(p_exp, hm.p));
            }
            { // hm_mt
                hist::ScatteringHistogram hm_mt = hist::HistogramManagerMT(&protein).calculate_all();
                REQUIRE(compare_hist(p_exp, hm_mt.p));
            }
            { // phm
                hist::ScatteringHistogram phm = protein.get_histogram();
                REQUIRE(compare_hist(p_exp, phm.p));
            }
            { // phm_mt
                hist::ScatteringHistogram phm_mt = hist::PartialHistogramManagerMT(&protein).calculate_all();
                REQUIRE(compare_hist(p_exp, phm_mt.p));
            }
        }

        SECTION("waters only") {
            // the following just describes the eight corners of a cube centered at origo, with an additional atom at the very middle
            vector<Atom> a = {};
            vector<Water> w = {Water(Vector3<double>(-1, -1, -1), 1, "C", "C", 1), Water(Vector3<double>(-1, 1, -1), 1, "C", "C", 1), 
                                Water(Vector3<double>(1, -1, -1), 1, "C", "C", 1),  Water(Vector3<double>(1, 1, -1), 1, "C", "C", 1), 
                                Water(Vector3<double>(-1, -1, 1), 1, "C", "C", 1),  Water(Vector3<double>(-1, 1, 1), 1, "C", "C", 1),
                                Water(Vector3<double>(1, -1, 1), 1, "C", "C", 1),   Water(Vector3<double>(1, 1, 1), 1, "C", "C", 1)};
            Protein protein(a, w);

            // set the weights to 1 so we can analytically determine the result
            for (auto& atom : protein.hydration_atoms) {
                atom.set_effective_charge(1);
            }
            protein.updated_charge = true;
            const vector<double> p_exp = {8, 0, 2*8*3, 8, 0, 0, 0, 0, 0, 0};

            { // hm
                hist::ScatteringHistogram hm = hist::HistogramManager(&protein).calculate_all();
                REQUIRE(compare_hist(p_exp, hm.p));
            }
            { // hm_mt
                hist::ScatteringHistogram hm_mt = hist::HistogramManagerMT(&protein).calculate_all();
                REQUIRE(compare_hist(p_exp, hm_mt.p));
            }
            { // phm
                hist::ScatteringHistogram phm = protein.get_histogram();
                REQUIRE(compare_hist(p_exp, phm.p));
            }
            { // phm_mt
                hist::ScatteringHistogram phm_mt = hist::PartialHistogramManagerMT(&protein).calculate_all();
                REQUIRE(compare_hist(p_exp, phm_mt.p));
            }
        }

        SECTION("both waters and atoms") {
            // the following just describes the eight corners of a cube centered at origo, with an additional atom at the very middle
            vector<Atom> b1 = {Atom(Vector3<double>(-1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, "C", "C", 1)};
            vector<Atom> b2 = {Atom(Vector3<double>(1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>(1, 1, -1), 1, "C", "C", 1)};
            vector<Atom> b3 = {Atom(Vector3<double>(-1, -1, 1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1, 1), 1, "C", "C", 1)};
            vector<Water> w = {Water(Vector3<double>(1, -1, 1), 1, "C", "C", 1),   Water(Vector3<double>(1, 1, 1), 1, "C", "C", 1)};
            vector<vector<Atom>> a = {b1, b2, b3};
            Protein protein(a, w);

            // set the weights to 1 so we can analytically determine the result
            // waters
            for (auto& atom : protein.hydration_atoms) {
                atom.set_effective_charge(1);
            }
            // atoms
            for (auto& body : protein.bodies) {
                for (auto& atom : body.atoms()) {
                    atom.set_effective_charge(1);
                }
            }
            protein.updated_charge = true;
            const vector<double> p_exp = {8, 0, 2*8*3, 8, 0, 0, 0, 0, 0, 0};

            { // hm
                hist::ScatteringHistogram hm = hist::HistogramManager(&protein).calculate_all();
                REQUIRE(compare_hist(p_exp, hm.p));
            }
            { // hm_mt
                hist::ScatteringHistogram hm_mt = hist::HistogramManagerMT(&protein).calculate_all();
                REQUIRE(compare_hist(p_exp, hm_mt.p));
            }
            { // phm
                hist::ScatteringHistogram phm = protein.get_histogram();
                REQUIRE(compare_hist(p_exp, phm.p));
            }
            { // phm_mt
                hist::ScatteringHistogram phm_mt = hist::PartialHistogramManagerMT(&protein).calculate_all();
                REQUIRE(compare_hist(p_exp, phm_mt.p));
            }
        }
    }

    SECTION("real data") {
        setting::protein::use_effective_charge = true;
                
        // create the atom, and perform a sanity check on our extracted list
        Protein protein("test/files/2epe.pdb");
        protein.generate_new_hydration();

        auto p_exp = protein.get_histogram().p;

        { // hm
            hist::ScatteringHistogram hm = hist::HistogramManager(&protein).calculate_all();
            REQUIRE(p_exp.size() == hm.p.size());
            REQUIRE(compare_hist(p_exp, hm.p));
        }
        { // hm_mt
            hist::ScatteringHistogram hm_mt = hist::HistogramManagerMT(&protein).calculate_all();
            REQUIRE(p_exp.size() == hm_mt.p.size());
            REQUIRE(compare_hist(p_exp, hm_mt.p));
        }
        { // phm
            hist::ScatteringHistogram phm = hist::PartialHistogramManager(&protein).calculate_all();
            REQUIRE(p_exp.size() == phm.p.size());
            REQUIRE(compare_hist(p_exp, phm.p));
        }
        { // phm_mt
            hist::ScatteringHistogram phm_mt = hist::PartialHistogramManagerMT(&protein).calculate_all();
            REQUIRE(p_exp.size() == phm_mt.p.size());
            REQUIRE(compare_hist(p_exp, phm_mt.p));
        }
    }
}

TEST_CASE("get_cm", "[protein]") {
    // make the protein
    vector<Atom> b1 = {Atom(1, "C", "", "LYS", "", 1, "", Vector3<double>(-1, -1, -1), 1, 0, "C", "0"),  Atom(2, "C", "", "LYS", "", 1, "", Vector3<double>(-1, 1, -1), 1, 0, "C", "0")};
    vector<Atom> b2 = {Atom(3, "C", "", "LYS", "", 1, "", Vector3<double>(1, -1, -1), 1, 0, "C", "0"), Atom(4, "C", "", "LYS", "", 1, "", Vector3<double>(1, 1, -1), 1, 0, "C", "0")};
    vector<Atom> b3 = {Atom(5, "C", "", "LYS", "", 1, "", Vector3<double>(-1, -1, 1), 1, 0, "C", "0"), Atom(6, "C", "", "LYS", "", 1, "", Vector3<double>(-1, 1, 1), 1, 0, "C", "0")};
    vector<Atom> b4 = {Atom(7, "C", "", "LYS", "", 1, "", Vector3<double>(1, -1, 1), 1, 0, "C", "0"),  Atom(8, "C", "", "LYS", "", 1, "", Vector3<double>(1, 1, 1), 1, 0, "C", "0")};
    vector<vector<Atom>> ap = {b1, b2, b3, b4};
    Protein protein(ap, {});

    Vector3<double> cm = protein.get_cm();
    REQUIRE(cm == Vector3<double>{0, 0, 0});
}

TEST_CASE("get_volume", "[protein]") {
    // make the protein
    vector<Atom> b1 = {Atom(1, "C", "", "LYS", "", 1, "", Vector3<double>(-1, -1, -1), 1, 0, "C", "0"),  Atom(2, "C", "", "LYS", "", 1, "", Vector3<double>(-1, 1, -1), 1, 0, "C", "0")};
    vector<Atom> b2 = {Atom(3, "C", "", "LYS", "", 1, "", Vector3<double>(1, -1, -1), 1, 0, "C", "0"), Atom(4, "C", "", "LYS", "", 1, "", Vector3<double>(1, 1, -1), 1, 0, "C", "0")};
    vector<Atom> b3 = {Atom(5, "C", "", "LYS", "", 1, "", Vector3<double>(-1, -1, 1), 1, 0, "C", "0"), Atom(6, "C", "", "LYS", "", 1, "", Vector3<double>(-1, 1, 1), 1, 0, "C", "0")};
    vector<Atom> b4 = {Atom(7, "C", "", "LYS", "", 1, "", Vector3<double>(1, -1, 1), 1, 0, "C", "0"),  Atom(8, "C", "", "LYS", "", 1, "", Vector3<double>(1, 1, 1), 1, 0, "C", "0")};
    vector<vector<Atom>> ap = {b1, b2, b3, b4};
    Protein protein(ap, {});

    REQUIRE_THAT(protein.get_volume_acids(), Catch::Matchers::WithinRel(4*constants::volume::amino_acids.get("LYS")));
}