#include "constants/ConstantsFwd.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <vector>
#include <string>
#include <iostream>

#include <data/Molecule.h>
#include <data/Body.h>
#include <data/record/Water.h>
#include <hydrate/Grid.h>
#include <hydrate/GridMember.h>
#include <constants/Constants.h>
#include <utility/Utility.h>
#include <fitter/LinearFitter.h>
#include <settings/All.h>
#include <fitter/HydrationFitter.h>
#include <hist/distance_calculator/HistogramManager.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>

using namespace data;
using namespace data::record;
using std::cout, std::endl, std::vector;

struct fixture {
    fixture() {
        settings::molecule::center = false;
        settings::molecule::use_effective_charge = false;
        settings::molecule::implicit_hydrogens = false;
    }

    Atom a1 = Atom(1, "C", "", "LYS", 'A', 1, "", Vector3<double>(-1, -1, -1), 1, 0, constants::atom_t::C, "0");
    Atom a2 = Atom(2, "C", "", "LYS", 'A', 1, "", Vector3<double>(-1,  1, -1), 1, 0, constants::atom_t::C, "0");
    Atom a3 = Atom(3, "C", "", "LYS", 'A', 1, "", Vector3<double>( 1, -1, -1), 1, 0, constants::atom_t::C, "0");
    Atom a4 = Atom(4, "C", "", "LYS", 'A', 1, "", Vector3<double>( 1,  1, -1), 1, 0, constants::atom_t::C, "0");
    Atom a5 = Atom(5, "C", "", "LYS", 'A', 1, "", Vector3<double>(-1, -1,  1), 1, 0, constants::atom_t::C, "0");
    Atom a6 = Atom(6, "C", "", "LYS", 'A', 1, "", Vector3<double>(-1,  1,  1), 1, 0, constants::atom_t::C, "0");
    Atom a7 = Atom(7, "C", "", "LYS", 'A', 1, "", Vector3<double>( 1, -1,  1), 1, 0, constants::atom_t::C, "0");
    Atom a8 = Atom(8, "C", "", "LYS", 'A', 1, "", Vector3<double>( 1,  1,  1), 1, 0, constants::atom_t::C, "0");
    
    Water w1 = Water(1, "O", "", "HOH", 'A', 1, "", Vector3<double>(-1, -1, -1), 1, 0, constants::atom_t::O, "0");
    Water w2 = Water(2, "O", "", "HOH", 'A', 1, "", Vector3<double>(-1,  1, -1), 1, 0, constants::atom_t::O, "0");    

    vector<Atom> b1 = {a1, a2};
    vector<Atom> b2 = {a3, a4};
    vector<Atom> b3 = {a5, a6};
    vector<Atom> b4 = {a7, a8};
    vector<Body> bodies = {Body(b1), Body(b2), Body(b3), Body(b4)};
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

TEST_CASE_METHOD(fixture, "Protein::Protein") {
    settings::molecule::use_effective_charge = false;
    settings::general::verbose = false;
    SECTION("Protein&") {
        Molecule protein1(bodies);
        Molecule protein2(protein1);
        REQUIRE(protein1.get_atoms() == protein2.get_atoms());
    }

    SECTION("vector<Body>&&") {
        Molecule protein(std::move(bodies));
        REQUIRE(protein.body_size() == 4);
        CHECK(protein.get_body(0).atom_size() == 2);
        CHECK(protein.get_body(1).atom_size() == 2);
        CHECK(protein.get_body(2).atom_size() == 2);
        CHECK(protein.get_body(3).atom_size() == 2);
    }

    SECTION("vector<string>&") {
        std::vector<std::string> files = {"test/files/2epe.pdb", "test/files/2epe.pdb"};
        Molecule protein(files);
        Body body("test/files/2epe.pdb"); // compare with body constructor
        body.get_waters().clear();

        REQUIRE(protein.body_size() == 2);
        CHECK(protein.atom_size() == 2002);
        CHECK(protein.water_size() == 96);
        CHECK(protein.get_body(0).equals_content(protein.get_body(1)));
        CHECK(protein.get_body(0).equals_content(body));
    }

    SECTION("ExistingFile&") {
        io::ExistingFile file("test/files/2epe.pdb");
        Molecule protein(file);
        Body body("test/files/2epe.pdb"); // compare with body constructor
        body.get_waters().clear();

        REQUIRE(protein.body_size() == 1);
        CHECK(protein.atom_size() == 1001);
        CHECK(protein.water_size() == 48);
        CHECK(protein.get_body(0).equals_content(body));
    }

    SECTION("vector<Body>&, vector<Water>&") {
        Molecule protein(bodies, {w1, w2});
        REQUIRE(protein.body_size() == 4);
        CHECK(protein.get_body(0).atom_size() == 2);
        CHECK(protein.get_body(1).atom_size() == 2);
        CHECK(protein.get_body(2).atom_size() == 2);
        CHECK(protein.get_body(3).atom_size() == 2);
        REQUIRE(protein.water_size() == 2);
        CHECK(protein.get_water(0) == w1);
        CHECK(protein.get_water(1) == w2);
    }

    SECTION("vector<Body>&") {
        Molecule protein(bodies);
        REQUIRE(protein.body_size() == 4);
        CHECK(protein.get_body(0).atom_size() == 2);
        CHECK(protein.get_body(1).atom_size() == 2);
        CHECK(protein.get_body(2).atom_size() == 2);
        CHECK(protein.get_body(3).atom_size() == 2);
    }

    SECTION("vector<Atom>&, vector<Water>&") {
        Molecule protein({a1, a2, a3, a4, a5, a6, a7, a8}, {w1, w2});
        REQUIRE(protein.body_size() == 1);
        CHECK(protein.get_body(0).atom_size() == 8);
        REQUIRE(protein.water_size() == 2);
        CHECK(protein.get_water(0) == w1);
        CHECK(protein.get_water(1) == w2);
    }

    SECTION("vector<Atom>&") {
        Molecule protein({a1, a2, a3, a4, a5, a6, a7, a8});
        REQUIRE(protein.body_size() == 1);
        CHECK(protein.get_body(0).atom_size() == 8);
    }
}

TEST_CASE("Protein::simulate_dataset") {
    settings::axes::qmax = 0.4;
    settings::molecule::use_effective_charge = false;
    settings::em::sample_frequency = 2;
    Molecule protein("test/files/2epe.pdb");

    SimpleDataset data = protein.simulate_dataset();
    fitter::LinearFitter fitter(data, protein.get_histogram());
    auto res = fitter.fit();
    REQUIRE_THAT(res->fval/res->dof, Catch::Matchers::WithinAbs(1., 0.5));
    // plots::PlotIntensityFit plot1(res);
    // plot1.save("figures/test/protein/check_chi2_1.png");
}

TEST_CASE_METHOD(fixture, "Protein::get_cm") {
    Molecule protein(bodies, {});
    Vector3<double> cm = protein.get_cm();
    REQUIRE(cm == Vector3<double>{0, 0, 0});
}

TEST_CASE_METHOD(fixture, "Protein::get_volume", "[broken]") {
    // broken since it was supposed to use the old protein.get_volume_acids() method
    // since the protein does not consist of a complete amino acid, the volume is not correct
    // TODO: create a protein containing a full amino acid and check if the volume is roughly correct
    Molecule protein(bodies, {});
    REQUIRE_THAT(protein.get_volume_grid(), Catch::Matchers::WithinRel(4*constants::volume::amino_acids.get("LYS")));
}

TEST_CASE_METHOD(fixture, "Protein::update_effective_charge") {
    settings::molecule::use_effective_charge = false;
    Molecule protein(bodies, {});

    double charge = protein.total_atomic_charge();
    double effective_charge = protein.total_effective_charge();
    REQUIRE(charge == effective_charge);

    protein.update_effective_charge(0.5);
    effective_charge = protein.total_effective_charge();
    REQUIRE(charge != effective_charge);

    protein.update_effective_charge(0);
    REQUIRE(charge == protein.total_effective_charge());
}

TEST_CASE_METHOD(fixture, "Protein::get_histogram") {
    SECTION("delegated to HistogramManager") {
        Molecule protein(bodies, {});
        REQUIRE(compare_hist(protein.get_histogram()->get_total_counts(), protein.get_histogram_manager()->calculate()->get_total_counts()));
    }
 
    SECTION("compare_debye") {
        vector<Atom> atoms = {Atom(Vector3<double>(-1, -1, -1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, constants::atom_t::C, "C", 1),
                              Atom(Vector3<double>( 1, -1, -1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>( 1, 1, -1), 1, constants::atom_t::C, "C", 1),
                              Atom(Vector3<double>(-1, -1,  1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>(-1, 1,  1), 1, constants::atom_t::C, "C", 1),
                              Atom(Vector3<double>( 1, -1,  1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>( 1, 1,  1), 1, constants::atom_t::C, "C", 1)};
        Molecule protein(atoms, {});

        vector<double> I_dumb = protein.debye_transform();
        vector<double> I_smart = protein.get_histogram()->debye_transform().get_counts();

        for (int i = 0; i < 8; i++) {
            if (!utility::approx(I_dumb[i], I_smart[i], 1e-1)) {
                cout << "Failed on index " << i << ". Values: " << I_dumb[i] << ", " << I_smart[i] << endl;
                REQUIRE(false);
            }
        }
        SUCCEED();
    }

    SECTION("compare_debye_real") {
        Molecule protein("test/files/2epe.pdb");
        protein.clear_hydration();

        std::cout << "hydration atoms: " << protein.get_waters().size() << std::endl; 

        vector<double> I_dumb = protein.debye_transform();
        vector<double> I_smart = protein.get_histogram()->debye_transform().get_counts();

        for (int i = 0; i < 8; i++) {
            if (!utility::approx(I_dumb[i], I_smart[i], 1e-3, 0.05)) {
                cout << "Failed on index " << i << ". Values: " << I_dumb[i] << ", " << I_smart[i] << endl;
                REQUIRE(false);
            }
        }
        SUCCEED();
    }
}

TEST_CASE_METHOD(fixture, "Protein::get_total_histogram") {
    Molecule protein(bodies, {});
    REQUIRE(compare_hist(protein.get_histogram()->get_total_counts(), protein.get_histogram_manager()->calculate_all()->get_total_counts()));
    // REQUIRE(protein.get_histogram() == protein.get_histogram_manager()->calculate_all());
}

TEST_CASE("Protein::save") {
    Molecule protein("test/files/2epe.pdb");
    protein.save("temp/test/protein_save_2epe.pdb");
    Molecule protein2("temp/test/protein_save_2epe.pdb");
    auto atoms1 = protein.get_atoms();
    auto atoms2 = protein2.get_atoms();
    REQUIRE(atoms1.size() == atoms2.size());
    for (size_t i = 0; i < atoms1.size(); ++i) {
        REQUIRE(atoms1[i].equals_content(atoms2[i]));
    }
}

TEST_CASE_METHOD(fixture, "Protein::generate_new_hydration") {
    settings::molecule::use_effective_charge = false;
    settings::general::verbose = false;

    SECTION("generates new waters") {
        // the molecule is really small, so we have to make sure there's enough space for the waters
        settings::grid::scaling = 5;
        Molecule protein(bodies);
        protein.generate_new_hydration();
        REQUIRE(protein.water_size() != 0);
        settings::grid::scaling = 0.25;
    }

    // we want to check that the hydration shells are consistent for fitting purposes
    SECTION("consistent hydration generation") {
        Molecule protein("test/files/2epe.pdb");
        fitter::LinearFitter fitter("test/files/2epe.dat", protein.get_histogram());

        protein.generate_new_hydration();
        double chi2 = fitter.fit()->fval;

        for (int i = 0; i < 10; i++) {
            protein.generate_new_hydration();
            double _chi2 = fitter.fit()->fval;
            REQUIRE_THAT(chi2, Catch::Matchers::WithinRel(_chi2));
        }
    }
}

TEST_CASE("Protein::get_volume_grid") {
    Molecule protein("test/files/2epe.pdb");
    REQUIRE(protein.get_volume_grid() == protein.get_grid()->get_volume());
}

// TEST_CASE("Protein::get_volume_calpha") {    
//     CHECK(false);
// }

TEST_CASE("Protein::molar_mass") {
    Molecule protein("test/files/2epe.pdb");
    REQUIRE(protein.molar_mass() == protein.absolute_mass()*constants::Avogadro);
}

TEST_CASE("Protein::absolute_mass") {
    Molecule protein("test/files/2epe.pdb");
    double sum = 0;
    for (auto& atom : protein.get_atoms()) {
        sum += atom.get_mass();
    }
    REQUIRE(protein.absolute_mass() == sum);
}

TEST_CASE("Protein::total_atomic_charge") {
    Molecule protein("test/files/2epe.pdb");
    double sum = 0;
    for (auto& atom : protein.get_atoms()) {
        sum += atom.get_absolute_charge();
    }
    REQUIRE(protein.total_atomic_charge() == sum);
}

TEST_CASE("Protein::total_effective_charge") {
    Molecule protein("test/files/2epe.pdb");
    double sum = 0;
    for (auto& atom : protein.get_atoms()) {
        sum += atom.get_effective_charge();
    }
    REQUIRE(protein.total_effective_charge() == sum);
}

TEST_CASE("Protein::get_relative_charge_density") {
    Molecule protein("test/files/2epe.pdb");
    REQUIRE_THAT(
        protein.get_relative_charge_density(), 
        Catch::Matchers::WithinAbs(
            (protein.total_atomic_charge() - constants::charge::density::water*protein.get_volume_grid())/protein.get_volume_grid(), 
            1e-6
        )
    );
}

TEST_CASE("Protein::get_relative_mass_density") {
    Molecule protein("test/files/2epe.pdb");
    REQUIRE(protein.get_relative_mass_density() == (protein.absolute_mass() - constants::mass::density::water*protein.get_volume_grid())/protein.get_volume_grid());
}

TEST_CASE("Protein::get_relative_charge") {
    Molecule protein("test/files/2epe.pdb");
    REQUIRE(protein.get_relative_charge() == protein.total_effective_charge() - protein.get_volume_grid()*constants::charge::density::water);
}

TEST_CASE_METHOD(fixture, "Protein::get_grid") {
    Molecule protein(bodies, {});
    // we just want to test that the grid is created by default
    REQUIRE(protein.get_grid() != nullptr);
}

TEST_CASE_METHOD(fixture, "Protein::set_grid") {
    Molecule protein(bodies, {});
    grid::Grid grid(Limit3D(0, 1, 0, 1, 0, 1));
    protein.set_grid(grid);
    REQUIRE(*protein.get_grid() == grid);
}

TEST_CASE_METHOD(fixture, "Protein::clear_hydration") {
    Molecule protein2(bodies, {w1, w2});
    REQUIRE(protein2.water_size() != 0);
    protein2.clear_hydration();
    REQUIRE(protein2.water_size() == 0);
}

TEST_CASE_METHOD(fixture, "Protein::center") {
    Molecule protein(bodies, {});
    REQUIRE(protein.get_cm() == Vector3<double>{0, 0, 0});

    protein.translate(Vector3<double>{1, 1, 1});
    REQUIRE(protein.get_cm() == Vector3<double>{1, 1, 1});
    
    protein.center();
    REQUIRE(protein.get_cm() == Vector3<double>{0, 0, 0});
}

TEST_CASE_METHOD(fixture, "Protein::get_body") {
    Molecule protein(bodies, {});
    REQUIRE(protein.get_body(0) == protein.get_bodies()[0]);
    REQUIRE(protein.get_body(1) == protein.get_bodies()[1]);
    REQUIRE(protein.get_body(2) == protein.get_bodies()[2]);
    REQUIRE(protein.get_body(3) == protein.get_bodies()[3]);
}

TEST_CASE_METHOD(fixture, "Protein::get_bodies") {
    Molecule protein(bodies, {});
    REQUIRE(protein.get_bodies() == bodies);
}

TEST_CASE_METHOD(fixture, "Protein::get_atoms") {
    Molecule protein(bodies, {});
    REQUIRE(protein.get_atoms() == vector<Atom>{a1, a2, a3, a4, a5, a6, a7, a8});
}

TEST_CASE_METHOD(fixture, "Protein::get_waters") {
    Molecule protein2(bodies, {w1, w2});
    REQUIRE(protein2.get_waters() == vector<Water>{w1, w2});
}

TEST_CASE_METHOD(fixture, "Protein::get_water") {
    Molecule protein2(bodies, {w1, w2});
    REQUIRE(protein2.get_water(0) == w1);
    REQUIRE(protein2.get_water(1) == w2);
}

TEST_CASE_METHOD(fixture, "Protein::create_grid") {
    Molecule protein(bodies, {});
    protein.clear_grid();
    auto grid = protein.get_grid();
    protein.create_grid();
    REQUIRE(protein.get_grid() != grid);
}

TEST_CASE_METHOD(fixture, "Protein::body_size") {
    Molecule protein(bodies, {});
    CHECK(protein.body_size() == 4);
}

TEST_CASE_METHOD(fixture, "Protein::atom_size") {
    Molecule protein(bodies, {});
    CHECK(protein.atom_size() == 8);
}

TEST_CASE_METHOD(fixture, "Protein::water_size") {
    Molecule protein(bodies, {});
    CHECK(protein.water_size() == 0);
    Molecule protein2(bodies, {w1, w2});
    CHECK(protein2.water_size() == 2);
}

TEST_CASE("Protein::fit") {
    Molecule protein("test/files/2epe.pdb");
    std::string measurement = "test/files/2epe.dat";
    fitter::HydrationFitter fitter(measurement, protein.get_histogram());

    auto pfit = protein.fit(measurement);
    auto hfit = fitter.fit();
    CHECK(pfit->fval == hfit->fval);
    CHECK(pfit->fevals == hfit->fevals);
    CHECK(pfit->parameters == hfit->parameters);
}

TEST_CASE_METHOD(fixture, "Protein::get_histogram_manager") {
    Molecule protein(bodies, {});
    CHECK(protein.get_histogram_manager() != nullptr);
}

// TEST_CASE_METHOD(fixture, "Protein::set_histogram_manager") {
//     Protein protein = Protein(bodies, {});
//     auto hm = protein.get_histogram_manager();
// }

TEST_CASE("Protein::translate") {
    Molecule protein("test/files/2epe.pdb");
    Vector3<double> cm = protein.get_cm();
    protein.translate(Vector3<double>{1, 1, 1});
    REQUIRE(protein.get_cm() == cm + Vector3<double>{1, 1, 1});
}

TEST_CASE("histogram") {
    settings::molecule::use_effective_charge = false;

    SECTION("multiple bodies, simple") {
        // make the protein
        vector<Atom> b1 = {Atom(Vector3<double>(-1, -1, -1), 1, constants::atom_t::C, "LYS", 1), Atom(Vector3<double>(-1, 1, -1), 1, constants::atom_t::C, "LYS", 1)};
        vector<Atom> b2 = {Atom(Vector3<double>( 1, -1, -1), 1, constants::atom_t::C, "LYS", 1), Atom(Vector3<double>( 1, 1, -1), 1, constants::atom_t::C, "LYS", 1)};
        vector<Atom> b3 = {Atom(Vector3<double>(-1, -1,  1), 1, constants::atom_t::C, "LYS", 1), Atom(Vector3<double>(-1, 1,  1), 1, constants::atom_t::C, "LYS", 1)};
        vector<Atom> b4 = {Atom(Vector3<double>( 1, -1,  1), 1, constants::atom_t::C, "LYS", 1), Atom(Vector3<double>( 1, 1,  1), 1, constants::atom_t::C, "LYS", 1)};
        vector<Body> ap = {Body(b1), Body(b2), Body(b3), Body(b4)};
        Molecule many(ap, {});

        // make the body
        vector<Atom> ab = {Atom(Vector3<double>(-1, -1, -1), 1, constants::atom_t::C, "LYS", 1), Atom(Vector3<double>(-1, 1, -1), 1, constants::atom_t::C, "LYS", 1),
                           Atom(Vector3<double>( 1, -1, -1), 1, constants::atom_t::C, "LYS", 1), Atom(Vector3<double>( 1, 1, -1), 1, constants::atom_t::C, "LYS", 1),
                           Atom(Vector3<double>(-1, -1,  1), 1, constants::atom_t::C, "LYS", 1), Atom(Vector3<double>(-1, 1,  1), 1, constants::atom_t::C, "LYS", 1),
                           Atom(Vector3<double>( 1, -1,  1), 1, constants::atom_t::C, "LYS", 1), Atom(Vector3<double>( 1, 1,  1), 1, constants::atom_t::C, "LYS", 1)};
        Molecule one(ab, {});

        // create some water molecules
        vector<Water> ws(10);
        for (size_t i = 0; i < ws.size(); i++) {
            ws[i] = Water::create_new_water(Vector3<double>(i, i, i));
        }

        many.get_waters() = ws;
        one.get_waters() = ws;

        // we now have a protein consisting of three bodies with the exact same contents as a single body.
        // the idea is now to compare the ScatteringHistogram output from their distance calculations, since it
        // is far easier to do for the single body. 
        auto d_m = many.get_histogram();
        auto d_o = one.get_histogram();

        // direct access to the histogram data (only p is defined)
        const vector<double>& p_m = d_m->get_total_counts();
        const vector<double>& p_o = d_o->get_total_counts();

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
        Body body("test/files/2epe.pdb");
        body.center();
        
        // We iterate through the protein data from the body, and split it into multiple pieces of size 100.  
        vector<Body> patoms; // vector containing the pieces we split it into
        vector<Atom> p_current(100); // vector containing the current piece
        unsigned int index = 0;      // current index in p_current
        for (unsigned int i = 0; i < body.get_atoms().size(); i++) {
            p_current[index] = body.get_atom(i);
            index++;
            if (index == 100) { // if index is 100, reset to 0
                patoms.emplace_back(p_current);
                index = 0;
            }
        }

        // add the final few atoms to our list
        if (index != 0) {
            p_current.resize(index);
            patoms.emplace_back(p_current);
        }

        // create the atom, and perform a sanity check on our extracted list
        Molecule protein(patoms, {});
        vector<Atom> protein_atoms = protein.get_atoms();
        vector<Atom> body_atoms = body.get_atoms();

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
        auto d_p = protein.get_histogram();
        auto d_b = hist::HistogramManager<false>(&protein).calculate_all();

        // direct access to the histogram data (only p is defined)
        const vector<double>& p = d_p->get_total_counts();
        const vector<double>& b_tot = d_b->get_total_counts();

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
        vector<Atom> atoms = {Atom(Vector3<double>(-1, -1, -1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, constants::atom_t::C, "C", 1),
                              Atom(Vector3<double>( 1, -1, -1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>( 1, 1, -1), 1, constants::atom_t::C, "C", 1),
                              Atom(Vector3<double>(-1, -1,  1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>(-1, 1,  1), 1, constants::atom_t::C, "C", 1),
                              Atom(Vector3<double>( 1, -1,  1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>( 1, 1,  1), 1, constants::atom_t::C, "C", 1)};

        // new auto-scaling approach
        Molecule protein1(atoms);
        grid::Grid grid1(atoms);
        protein1.set_grid(grid1);

        // old approach
        Molecule protein2(atoms);
        grid::Grid grid2({-2, 2, -2, 2, -2, 2}); 
        grid2.add(atoms);
        protein2.set_grid(grid2);

        // generate the distance histograms
        auto h1 = protein1.get_histogram();
        auto h2 = protein2.get_histogram();

        // direct access to the histogram data (only p is defined)
        const vector<double>& p1 = h1->get_total_counts();
        const vector<double>& p2 = h2->get_total_counts();

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

// struct fixture {
//     Atom a1 = Atom(Vector3<double>(-1, -1, -1), 1, "C", "C", 1);
//     Atom a2 = Atom(Vector3<double>(-1,  1, -1), 1, "C", "C", 1);
//     Atom a3 = Atom(Vector3<double>(-1, -1,  1), 1, "C", "C", 1);
//     Atom a4 = Atom(Vector3<double>(-1,  1,  1), 1, "C", "C", 1);
//     Atom a5 = Atom(Vector3<double>( 1, -1, -1), 1, "C", "C", 1);
//     Atom a6 = Atom(Vector3<double>( 1,  1, -1), 1, "C", "C", 1);
//     Atom a7 = Atom(Vector3<double>( 1, -1,  1), 1, "C", "C", 1);
//     Atom a8 = Atom(Vector3<double>( 1,  1,  1), 1, "He", "He", 1);

//     Body b1 = Body(std::vector<Atom>{a1, a2});
//     Body b2 = Body(std::vector<Atom>{a3, a4});
//     Body b3 = Body(std::vector<Atom>{a5, a6});
//     Body b4 = Body(std::vector<Atom>{a7, a8});
//     std::vector<Body> ap = {b1, b2, b3, b4};
//     Protein protein = Protein(ap);
// };

#include <data/state/StateManager.h>
#include <data/state/BoundSignaller.h>
#include <hist/distance_calculator/HistogramManagerFactory.h>
TEST_CASE_METHOD(fixture, "Protein::bind_body_signallers") {
    Molecule protein(bodies, {});
    settings::general::verbose = false;

    SECTION("at construction") {
        auto& bodies = protein.get_bodies();
        REQUIRE(bodies.size() == 4);
        auto manager = protein.get_histogram_manager()->get_state_manager();
        for (unsigned int i = 0; i < bodies.size(); ++i) {
            CHECK(std::dynamic_pointer_cast<signaller::BoundSignaller>(bodies[i].get_signaller()) != nullptr);
            CHECK(manager->get_probe(i) == bodies[i].get_signaller());
        }

        manager->reset_to_false();
        for (unsigned int i = 0; i < bodies.size(); ++i) {
            bodies[i].changed_external_state();
            CHECK(manager->is_externally_modified(i));
        }
    }

    SECTION("after construction") {
        auto& bodies = protein.get_bodies();
        REQUIRE(bodies.size() == 4);
        protein.set_histogram_manager(hist::factory::construct_histogram_manager(&protein));
        auto manager = protein.get_histogram_manager()->get_state_manager();

        for (unsigned int i = 0; i < bodies.size(); ++i) {
            CHECK(std::dynamic_pointer_cast<signaller::BoundSignaller>(bodies[i].get_signaller()) != nullptr);
            CHECK(manager->get_probe(i) == bodies[i].get_signaller());
        }
    }
}

TEST_CASE_METHOD(fixture, "Protein::signal_modified_hydration_layer") {
    Molecule protein(bodies, {});
    auto manager = protein.get_histogram_manager()->get_state_manager();
    manager->reset_to_false();
    REQUIRE(manager->get_modified_hydration() == false);

    protein.signal_modified_hydration_layer();
    REQUIRE(manager->get_modified_hydration() == true);
}