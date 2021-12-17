// includes
#include <vector>
#include <string>
#include <fstream>

#include "tests/Test.h"
#include "data/Protein.h"
#include "hydrate/Grid.h"
#include "constants.h"

void create_test_file() {
    std::ofstream file1("temp1.pdb");
    std::ofstream file2("temp2.pdb");
    std::ofstream file3("temp3.pdb");
    // the following just describes the eight corners of a cube centered at origo, with an additional atom at the very middle
    file1 << "ATOM      1  C   LYS A   1          -1      -1      -1  1.00 00.00           C \n"
          << "ATOM      2  C   LYS A   1          -1       1      -1  1.00 00.00           C \n"
          << "ATOM      3  C   LYS A   1           1      -1      -1  1.00 00.00           C \n"
          << "ATOM      4  C   LYS A   1           1       1      -1  1.00 00.00           C \n";
    file1.close();
    file2 << "ATOM      5  C   LYS A   1          -1      -1       1  1.00 00.00           C \n"
          << "ATOM      6  C   LYS A   1          -1       1       1  1.00 00.00           C \n"
          << "ATOM      7  C   LYS A   1           1      -1       1  1.00 00.00           C \n"
          << "ATOM      8  C   LYS A   1           1       1       1  1.00 00.00           C \n";
    file2.close();
    file3 << "ATOM      9  C   LYS A   1           0       0       0  1.00 00.00           C  ";
    file3.close();
}

void test_calc_distances() {
    std::ofstream file("temp_dist.pdb");
    // the following just describes the eight corners of a cube centered at origo, with an additional atom at the very middle
    file << "ATOM      1  C   LYS A   1          -1      -1      -1  1.00 00.00           C \n"
         << "ATOM      2  C   LYS A   1          -1       1      -1  1.00 00.00           C \n"
         << "ATOM      3  C   LYS A   1           1      -1      -1  1.00 00.00           C \n"
         << "ATOM      4  C   LYS A   1           1       1      -1  1.00 00.00           C \n"

         << "ATOM      5  C   LYS A   1          -1      -1       1  1.00 00.00           C \n"
         << "ATOM      6  C   LYS A   1          -1       1       1  1.00 00.00           C \n"
         << "ATOM      7  C   LYS A   1           1      -1       1  1.00 00.00           C \n"
         << "ATOM      8  C   LYS A   1           1       1       1  1.00 00.00           C \n"

         << "ATOM      9  C   LYS A   1           0       0       0  1.00 00.00           C  ";
    file.close();

    // create some water molecules
    vector<Hetatom> atoms(10);
    for (size_t i = 0; i < atoms.size(); i++) {
        atoms[i] = Hetatom::create_new_water(Vector3(i, i, i));
    }

    vector<string> files = {"temp1.pdb", "temp2.pdb", "temp3.pdb"};
    Protein protein(files);
    // Protein protein("temp1.pdb");
    Body body("temp_dist.pdb");
    remove("temp_dist.pdb");

    body.hydration_atoms = atoms;
    protein.hydration_atoms = atoms;

    // we now have a protein consisting of three bodies with the exact same contents as a single body.
    // the idea is now to compare the ScatteringHistogram output from their distance calculations, since it
    // is far easier to do for the single body. 
    shared_ptr<ScatteringHistogram> d_b = body.get_histogram();
    shared_ptr<ScatteringHistogram> d_p = protein.get_histogram();

    // direct access to the histogram data (only p_tot is defined)
    const vector<double>& p_tot = d_p->p_tot;
    const vector<double>& b_tot = d_b->p_tot;

    // compare each entry
    for (size_t i = 0; i < b_tot.size(); i++) {
        if (!approx(p_tot[i], b_tot[i])) {
            IS_TRUE(false);
            cout << "Failed on index " << i << ". Values: " << p_tot[i] << ", " << b_tot[i] << endl;
            // break;
        }
    }
}

void test_calc_distances_2epe() {
    Body body("data/2epe.pdb");

    // We iterate through the protein data from the body, and split it into multiple pieces of size 100.  
    vector<vector<Atom>> patoms; // vector containing the pieces we split it into
    vector<Atom> p_current(100); // vector containing the current piece
    size_t index = 0; // current index in p_current
    for (size_t i = 0; i < body.protein_atoms.size(); i++) {
        p_current[index] = body.protein_atoms[i];
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
    vector<Atom> protein_atoms = protein.get_protein_atoms();
    vector<Atom> body_atoms = body.get_protein_atoms();

    // sizes must be equal. this also serves as a separate consistency check on the body generation. 
    if (protein_atoms.size() != body_atoms.size()) {
        IS_TRUE(false);
        cout << "Sizes " << protein_atoms.size() << " and " << body_atoms.size() << " should be equal. " << endl;
        return;
    }

    // stronger consistency check - we check that all atoms are equal, and appear in the exact same order
    for (size_t i = 0; i < protein_atoms.size(); i++) {
        if (protein_atoms[i] != body_atoms[i]) {
            IS_TRUE(false);
            cout << "Comparison failed on index " << i << endl;
            return;
        }
    }

    // generate a hydration layer for the body, and copy it over to the protein
    body.generate_new_hydration();
    protein.hydration_atoms = body.hydration_atoms;

    // generate the distance histograms
    shared_ptr<ScatteringHistogram> d_b = body.get_histogram();
    shared_ptr<ScatteringHistogram> d_p = protein.get_histogram();

    cout << "GRID SIZES: protein: " << protein.get_grid()->get_protein_atoms().size() << ", body: " << body.get_grid()->get_protein_atoms().size() << endl;

    // direct access to the histogram data (only p_tot is defined)
    const vector<double>& p_tot = d_p->p_tot;
    const vector<double>& b_tot = d_b->p_tot;

    // compare each entry
    for (size_t i = 0; i < b_tot.size(); i++) {
        if (!approx(p_tot[i], b_tot[i])) {
            IS_TRUE(false);
            cout << "Failed on index " << i << ". Values: " << p_tot[i] << ", " << b_tot[i] << endl;
            break;
        }
    }
}

void test_get_cm() {
    Protein protein({"temp1.pdb", "temp2.pdb", "temp3.pdb"});
    Vector3 cm = protein.get_cm();
    IS_TRUE(cm == Vector3({0, 0, 0}));
}

void test_volume() {
    Protein protein({"temp1.pdb", "temp2.pdb", "temp3.pdb"});
    IS_TRUE(protein.get_volume_acids() == 3*constants::volume::lysine);
}

int main(void) {
    setting::grid::psc = setting::grid::RadialStrategy;
    cout << "Summary of Protein testing:" << endl;
    create_test_file();
    // test_get_cm();
    // test_volume();
    // test_calc_distances();
    test_calc_distances_2epe();
    remove("temp1.pdb");
    remove("temp2.pdb");
    remove("temp3.pdb");

    if (passed_all) {
        cout << "\033[1;32m" << "All Protein tests passed." << "\033[0m" << endl;
    } else {
        print_err("Some Protein tests failed.");
    }
    return 0;
}