// includes
#include <vector>
#include <string>
#include <fstream>

// ROOT
#include <TStyle.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TH1.h>

// my own includes
#include "../Protein.cpp"

using namespace ROOT;

int main(int argc, char const *argv[])
{
    cout << "Testing Protein class...\t\r" << std::flush;
    std::ofstream file("temp.pdb");
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

    Protein protein("temp.pdb");

    // check get_cm
    TVector3 cm = protein.get_cm();
    try {
        assert(cm[0] == 0);
        assert(cm[1] == 0);
        assert(cm[2] == 0);
    } catch (const std::exception& e) {
        print_err("get_cm failed.");
    }

    // check the grid generator
    double width = 0.1;
    auto[corner, bins] = protein.generate_grid(width);
    try {
        assert(corner[0] == -1);
        assert(corner[1] == -1);
        assert(corner[2] == -1);
        assert(bins[0] == 2/width);
    } catch (const std::exception& e) {
        print_err("generate_grid failed.");
    }

    // check find_protein_locations
    // std::map<int, vector<int>> locs = protein.find_protein_locations();
    // vector<vector<vector<bool>>> grid = protein.boolean_representation(locs);
    // try { // these are just all corner locations
    //     assert(locs[0][0][0]);
    //     assert(locs[0][0][bins[2]]);
    //     assert(locs[0][bins[1]][0]);
    //     assert(locs[0][bins[1]][bins[2]]);
    //     assert(locs[bins[0]][0][0]);
    //     assert(locs[bins[0]][0][bins[2]]);
    //     assert(locs[bins[0]][bins[1]][0]);
    //     assert(locs[bins[0]][bins[1]][bins[2]]);
    // } catch (const std::exception& e) {
    //     print_err("find_protein_locations failed.");
    // }

    remove("temp.pdb");

    cout << "\033[1;32m" << "All protein tests passed." << "\033[0m" << endl;
    return 0;
}