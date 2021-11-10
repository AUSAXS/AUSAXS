// includes
#include <vector>
#include <string>
#include <fstream>

#include "tests/Test.h"
#include "Protein.cpp"
#include "Grid.cpp"

using namespace ROOT;

/**
 * @brief Test that the grid is generated correctly.
 */
void test_grid_generation() {
    TVector3 base(-10, -10, -10);
    int width = 1;
    int bins = 21;
    Grid grid(base, width, bins);

    shared_ptr<Atom> atom = std::make_shared<Atom>(Atom({0, 0, 0}, 0, "C", "", 0));
    vector<shared_ptr<Atom>> a = {atom};
    grid.add(&a);
    vector<vector<vector<bool>>> &g = *grid.get_grid();

    // check that it was placed correctly in the grid
    IS_TRUE(g[10][10][10]);
    IS_TRUE(!g[10][10][11]);
    IS_TRUE(!g[10][11][11]);
    IS_TRUE(!g[11][10][10]);
    IS_TRUE(!g[9][8][14]);
}

void test_simple_bounding_box() {
    TVector3 base(-10, -10, -10);
    int width = 1;
    int bins = 21;
    Grid grid(base, width, bins);

    shared_ptr<Atom> atom = std::make_shared<Atom>(Atom({0, 0, 0}, 0, "C", "", 0));
    vector<shared_ptr<Atom>> a = {atom};
    grid.add(&a);
    vector<vector<vector<bool>>> &g = *grid.get_grid();

    vector<vector<int>> box = grid.bounding_box();
    IS_TRUE(box[0][0] == 10);
    IS_TRUE(box[0][1] == 10);
    IS_TRUE(box[1][0] == 10);
    IS_TRUE(box[1][1] == 10);
    IS_TRUE(box[2][0] == 10);
    IS_TRUE(box[2][1] == 10);
}

void test_complex_bounding_box() {
    TVector3 base(-10, -10, -10);
    int width = 1;
    int bins = 21;
    Grid grid(base, width, bins, 3);

    shared_ptr<Atom> a1 = std::make_shared<Atom>(Atom({5, 0, -7}, 0, "C", "", 1));
    shared_ptr<Atom> a2 = std::make_shared<Atom>(Atom({0, -5, 0}, 0, "C", "", 2));
    vector<shared_ptr<Atom>> a = {a1, a2};

    grid.add(&a);
    grid.expand_volume();
    vector<vector<vector<bool>>> &g = *grid.get_grid();

    vector<vector<int>> box = grid.bounding_box();
    IS_TRUE(box[0][0] == 10);
    IS_TRUE(box[0][1] == 15);
    IS_TRUE(box[1][0] == 5);
    IS_TRUE(box[1][1] == 10);
    IS_TRUE(box[2][0] == 3);
    IS_TRUE(box[2][1] == 10);
}

void test_find_free_locs() {
    TVector3 base(-10, -10, -10);
    int width = 1;
    int bins = 21;
    Grid grid(base, width, bins, 3);

    shared_ptr<Atom> atom = std::make_shared<Atom>(Atom({0, 0, 0}, 0, "C", "", 0));
    vector<shared_ptr<Atom>> a = {atom};
    grid.add(&a);
    grid.expand_volume();
    vector<vector<vector<bool>>> &g = *grid.get_grid();

    vector<vector<int>> locs = grid.find_free_locs();
    IS_TRUE(locs.size() == 6);
    if (locs.size() == 6) { // avoid crashing if the above fails
        IS_TRUE(locs[0] == vector({4, 10, 10})); // (-2r, 0, 0)
        IS_TRUE(locs[1] == vector({10, 4, 10})); // (0, -2r, 0)
        IS_TRUE(locs[2] == vector({10, 10, 4})); // (0, 0, -2r)
        IS_TRUE(locs[3] == vector({10, 10, 16})); // (0, 0, 2r)
        IS_TRUE(locs[4] == vector({10, 16, 10})); // (0, 2r, 0)
        IS_TRUE(locs[5] == vector({16, 10, 10})); // (2r, 0, 0)
    }
}

void test_volume_expansion() {
    TVector3 base(-10, -10, -10);
    int width = 1;
    int bins = 21;
    int radius = 3;
    Grid grid(base, width, bins, radius);

    shared_ptr<Atom> atom = std::make_shared<Atom>(Atom({0, 0, 0}, 0, "C", "", 0));
    vector<shared_ptr<Atom>> a = {atom};
    grid.add(&a);
    vector<vector<vector<bool>>> &g = *grid.get_grid();

    grid.expand_volume();
    // check that the x=13 plane looks like this: 
    // x x x
    // x o x
    // x x x
    IS_FALSE(g[14][10][10]); // check that it does not extend to the x=14 plane
    IS_TRUE(g[13][10][10]);

    IS_FALSE(g[13][9][9]);
    IS_FALSE(g[13][9][10]);
    IS_FALSE(g[13][9][11]);

    IS_FALSE(g[13][10][9]);
    IS_FALSE(g[13][11][11]);

    IS_FALSE(g[13][11][9]);
    IS_FALSE(g[13][11][10]);
    IS_FALSE(g[13][11][11]);

    // repeat with the x=7 plane
    IS_FALSE(g[6][10][10]); // check that it does not extend to the x=6 plane
    IS_TRUE(g[7][10][10]);

    IS_FALSE(g[7][9][9]);
    IS_FALSE(g[7][9][10]);
    IS_FALSE(g[7][9][11]);

    IS_FALSE(g[7][10][9]);
    IS_FALSE(g[7][11][11]);

    IS_FALSE(g[7][11][9]);
    IS_FALSE(g[7][11][10]);
    IS_FALSE(g[7][11][11]);

    // check some other points as well
    // x=10, z=10 line, from z=6 it looks like x o o o o o o o x
    IS_FALSE(g[10][6][10]);
    IS_TRUE(g[10][7][10]);
    IS_TRUE(g[10][8][10]);
    IS_TRUE(g[10][9][10]);
    IS_TRUE(g[10][10][10]);
    IS_TRUE(g[10][11][10]);
    IS_TRUE(g[10][12][10]);
    IS_TRUE(g[10][13][10]);
    IS_FALSE(g[10][14][10]);

    // x=10, y=10 line, from y=6 it looks like x o o o o o o o x
    IS_FALSE(g[10][10][6]);
    IS_TRUE(g[10][10][7]);
    IS_TRUE(g[10][10][8]);
    IS_TRUE(g[10][10][9]);
    IS_TRUE(g[10][10][10]);
    IS_TRUE(g[10][10][11]);
    IS_TRUE(g[10][10][12]);
    IS_TRUE(g[10][10][13]);
    IS_FALSE(g[10][10][14]);

    // some random points
    IS_TRUE(g[9][9][9]);
    IS_TRUE(!g[8][8][8]);
    IS_TRUE(!g[13][13][13]);
}

int main(void)
{
    cout << "Summary of Grid testing:" << std::endl;
    test_grid_generation();
    test_simple_bounding_box();
    test_complex_bounding_box();
    test_find_free_locs();
    test_volume_expansion();

    if (passed_all) {
        cout << "\033[1;32m" << "All Grid tests passed." << "\033[0m" << endl;
    } else {
        print_err("Some Grid tests failed.");
    }
    return 0;
}