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
#include "../Grid.cpp"

using namespace ROOT;

int main(int argc, char const *argv[])
{
    cout << "Testing Grid class...\t\r" << std::flush;
    TVector3 base(-10, -10, -10);
    int width = 1;
    int bins = 21;
    Grid grid(base, width, bins);

    // check the grid generation
    Atom* atom = new Atom({0, 0, 0}, 0, "C", "");
    vector<Atom*> a = {atom};
    grid.add(&a);
    vector<vector<vector<bool>>> &g = *grid.get_grid();

    try {
        // check that it was placed correctly in the grid
        assert(g[10][10][10]);

        // some additional random checks
        assert(!g[10][10][11]);
        assert(!g[10][11][11]);
        assert(!g[11][10][10]);
        assert(!g[9][8][14]);
    } catch (const std::exception& e) {
        print_err("Atom was placed incorrectly in the Grid.");
    }

    // check that the bounding box of a point atom is correct
    vector<vector<int>> box = grid.bounding_box();
    try {
        // check that it was placed correctly in the grid
        assert(box[0][0] == 10);
        assert(box[0][1] == 10);
        assert(box[1][0] == 10);
        assert(box[1][1] == 10);
        assert(box[2][0] == 10);
        assert(box[2][1] == 10);
    } catch (const std::exception& e) {
        print_err("The bounding box is incorrect.");
    }

    // check the possible HOH spot locator
    grid.set_radius(3);
    vector<vector<int>> locs = grid.find_free_locs();
    try {
        assert(locs.size() == 6);
        assert(locs[0] == vector({7, 10, 10}));
        assert(locs[1] == vector({13, 10, 10}));
        assert(locs[2] == vector({10, 7, 10}));
        assert(locs[3] == vector({10, 13, 10}));
        assert(locs[4] == vector({10, 10, 7}));
        assert(locs[5] == vector({10, 10, 13}));
    } catch (const std::exception& e) {
        print_err("Volume expansion failed.");
    }

    // check that volumes are filled correctly
    grid.expand_volume(3);
    try {
        assert(g[10][10][12]);
        assert(g[12][10][10]);
        assert(g[10][12][10]);
        assert(g[9][9][9]);

        assert(!g[8][8][8]);
        assert(!g[12][12][12]);
    } catch (const std::exception& e) {
        print_err("Volume expansion failed.");
    }

    // check that the bounding box of a more advanced structure is correct
    Atom* a1 = new Atom({5, 0, -7}, 0, "C", "");
    Atom* a2 = new Atom({0, -5, 0}, 0, "C", "");
    a1->set_serial(1);
    a2->set_serial(2);
    a = {a1, a2};
    grid.add(&a);
    grid.expand_volume(3);
    box = grid.bounding_box();
    try {
        assert(box[0][0] == 10);
        assert(box[0][1] == 15);
        assert(box[1][0] == 5);
        assert(box[1][1] == 10);
        assert(box[2][0] == 3);
        assert(box[2][1] == 10);
    } catch (const std::exception& e) {
        print_err("The bounding box for three atoms is incorrect.");
    }

    cout << "\033[1;32m" << "All Grid tests passed." << "\033[0m" << endl;
    return 0;
}