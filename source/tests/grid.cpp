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
    Grid grid(base, width, bins);

    shared_ptr<Atom> a1 = std::make_shared<Atom>(Atom({5, 0, -7}, 0, "C", "", 1));
    shared_ptr<Atom> a2 = std::make_shared<Atom>(Atom({0, -5, 0}, 0, "C", "", 2));
    vector<shared_ptr<Atom>> a = {a1, a2};

    grid.add(&a);
    grid.expand_volume(3);
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
    Grid grid(base, width, bins);

    shared_ptr<Atom> atom = std::make_shared<Atom>(Atom({0, 0, 0}, 0, "C", "", 0));
    vector<shared_ptr<Atom>> a = {atom};
    grid.add(&a);
    vector<vector<vector<bool>>> &g = *grid.get_grid();

    grid.set_radius(3);
    vector<vector<int>> locs = grid.find_free_locs();
    IS_TRUE(locs.size() == 6);
    IS_TRUE(locs[0] == vector({7, 10, 10}));
    IS_TRUE(locs[1] == vector({13, 10, 10}));
    IS_TRUE(locs[2] == vector({10, 7, 10}));
    IS_TRUE(locs[3] == vector({10, 13, 10}));
    IS_TRUE(locs[4] == vector({10, 10, 7}));
    IS_TRUE(locs[5] == vector({10, 10, 13}));
}

void test_volume_expansion() {
    TVector3 base(-10, -10, -10);
    int width = 1;
    int bins = 21;
    Grid grid(base, width, bins);

    shared_ptr<Atom> atom = std::make_shared<Atom>(Atom({0, 0, 0}, 0, "C", "", 0));
    vector<shared_ptr<Atom>> a = {atom};
    grid.add(&a);
    vector<vector<vector<bool>>> &g = *grid.get_grid();

    grid.expand_volume(3);
    IS_TRUE(g[10][10][12]);
    IS_TRUE(g[12][10][10]);
    IS_TRUE(g[10][12][10]);
    IS_TRUE(g[9][9][9]);
    IS_TRUE(!g[8][8][8]);
    IS_TRUE(!g[12][12][12]);
}

int main(void)
{
    cout << "Summary of Grid testing:" << std::endl;
    test_grid_generation();
    test_simple_bounding_box();
    test_complex_bounding_box();
    test_find_free_locs();
    test_volume_expansion();
    cout << "\033[1;32m" << "All Grid tests passed." << "\033[0m" << endl;
    return 0;
}