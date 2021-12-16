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
    // file1 << "ATOM      1  C   LYS A   1          -1      -1      -1  1.00 00.00           C \n";
    // file1.close();
    // file2 << "ATOM      2  C   LYS A   1          -1       1      -1  1.00 00.00           C \n";
    // file2.close();
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
    // file << "ATOM      1  C   LYS A   1          -1      -1      -1  1.00 00.00           C \n";
    // file.close();

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

    // body.protein_atoms = {};

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

void test_calc_distances2() {
    // check if file was succesfully opened
    std::ifstream input("data/2epe.pdb");
    if (!input.is_open()) {return;}

    string line; // placeholder for the current line
    size_t counter = 1;
    vector<string> files;
    string out = "";
    while(getline(input, line)) {
        if (!(line.substr(0, 4) == "ATOM" || line.substr(0, 6) == "HETATOM")) {continue;}
        out += line;
        if (counter % 100 == 0) {
            string file_name = "temp" + std::to_string(files.size()+1) + ".pdb";
            files.push_back(file_name); 
            std::ofstream file(file_name);
            file << out;
            file.close();
            out = "";
        }
        counter++;
    }

    Protein protein(files);
    Body body("data/2epe.pdb");

    body.generate_new_hydration();
    protein.hydration_atoms = body.hydration_atoms;

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
    cout << "Summary of Protein testing:" << endl;
    create_test_file();
    // test_get_cm();
    // test_volume();
    test_calc_distances();
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