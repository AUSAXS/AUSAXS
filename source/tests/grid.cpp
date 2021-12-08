// includes
#include <vector>
#include <string>
#include <fstream>

#include "tests/Test.h"
#include "Protein.h"
#include "hydrate/Grid.h"
#include "settings.h"

#include <math/Vector3.h>
 
using namespace ROOT;

/**
 * @brief Test that the grid is generated correctly.
 */
void test_grid_generation() {
    Vector3 base(-10, -10, -10);
    int width = 1;
    int bins = 21;
    Grid grid(base, width, bins);

    vector<Atom> a = {Atom({0, 0, 0}, 0, "C", "", 0)};
    grid.add(a);
    vector<vector<vector<char>>> &g = grid.grid;

    // check that it was placed correctly in the grid
    IS_TRUE(g[10][10][10] == 'A');
    IS_TRUE(g[10][10][11] == 0);
    IS_TRUE(g[10][11][11] == 0);
    IS_TRUE(g[11][10][10] == 0);
    IS_TRUE(g[9][8][14] == 0);
}

void test_simple_bounding_box() {
    Vector3 base(-10, -10, -10);
    int width = 1;
    int bins = 21;
    Grid grid(base, width, bins);

    vector<Atom> a = {Atom({0, 0, 0}, 0, "C", "", 0)};
    grid.add(a);

    vector<vector<int>> box = grid.bounding_box();
    IS_TRUE(box[0][0] == 10);
    IS_TRUE(box[0][1] == 10);
    IS_TRUE(box[1][0] == 10);
    IS_TRUE(box[1][1] == 10);
    IS_TRUE(box[2][0] == 10);
    IS_TRUE(box[2][1] == 10);
}

void test_complex_bounding_box() {
    Vector3 base(-10, -10, -10);
    int width = 1;
    int bins = 21;
    Grid grid(base, width, bins, 3);

    vector<Atom> a = {Atom({5, 0, -7}, 0, "C", "", 1), Atom({0, -5, 0}, 0, "C", "", 2)};
    grid.add(a);
    grid.expand_volume();

    vector<vector<int>> box = grid.bounding_box();
    IS_TRUE(box[0][0] == 10);
    IS_TRUE(box[0][1] == 15);
    IS_TRUE(box[1][0] == 5);
    IS_TRUE(box[1][1] == 10);
    IS_TRUE(box[2][0] == 3);
    IS_TRUE(box[2][1] == 10);
}

void test_volume_expansion() {
    Vector3 base(-10, -10, -10);
    int width = 1;
    int bins = 21;
    int radius = 3;
    Grid grid(base, width, bins, radius);

    vector<Atom> a = {Atom({0, 0, 0}, 0, "C", "", 0)};
    grid.add(a);
    vector<vector<vector<char>>> &g = grid.grid;
    grid.expand_volume();

    IS_TRUE(g[10][10][10] == 'A'); // check that the center is still marked as capital 'A'

    // check that the x=13 plane looks like this: 
    // x x x
    // x o x
    // x x x
    IS_TRUE(g[14][10][10] == 0); // check that it does not extend to the x=14 plane
    IS_TRUE(g[13][10][10] == 'a');

    IS_TRUE(g[13][9][9] == 0);
    IS_TRUE(g[13][9][10] == 0);
    IS_TRUE(g[13][9][11] == 0);

    IS_TRUE(g[13][10][9] == 0);
    IS_TRUE(g[13][11][11] == 0);

    IS_TRUE(g[13][11][9] == 0);
    IS_TRUE(g[13][11][10] == 0);
    IS_TRUE(g[13][11][11] == 0);

    // repeat with the x=7 plane
    IS_TRUE(g[6][10][10] == 0); // check that it does not extend to the x=6 plane
    IS_TRUE(g[7][10][10] == 'a');

    IS_TRUE(g[7][9][9] == 0);
    IS_TRUE(g[7][9][10] == 0);
    IS_TRUE(g[7][9][11] == 0);

    IS_TRUE(g[7][10][9] == 0);
    IS_TRUE(g[7][11][11] == 0);

    IS_TRUE(g[7][11][9] == 0);
    IS_TRUE(g[7][11][10] == 0);
    IS_TRUE(g[7][11][11] == 0);

    // check some other points as well
    // x=10, z=10 line, from z=6 it looks like x o o o o o o o x
    IS_TRUE(g[10][6][10] == 0);
    IS_TRUE(g[10][7][10] == 'a');
    IS_TRUE(g[10][8][10] == 'a');
    IS_TRUE(g[10][9][10] == 'a');
    IS_TRUE(g[10][10][10] == 'A');
    IS_TRUE(g[10][11][10] == 'a');
    IS_TRUE(g[10][12][10] == 'a');
    IS_TRUE(g[10][13][10] == 'a');
    IS_TRUE(g[10][14][10] == 0);

    // x=10, y=10 line, from y=6 it looks like x o o o o o o o x
    IS_TRUE(g[10][10][6] == 0);
    IS_TRUE(g[10][10][7] == 'a');
    IS_TRUE(g[10][10][8] == 'a');
    IS_TRUE(g[10][10][9] == 'a');
    IS_TRUE(g[10][10][10] == 'A');
    IS_TRUE(g[10][10][11] == 'a');
    IS_TRUE(g[10][10][12] == 'a');
    IS_TRUE(g[10][10][13] == 'a');
    IS_TRUE(g[10][10][14] == 0);

    // some random points
    IS_TRUE(g[9][9][9] == 'a');
    IS_TRUE(g[8][8][8] == 0);
    IS_TRUE(g[13][13][13] == 0);
}

void test_hydrate() {
    Vector3 base(-10, -10, -10);
    double width = 1;
    int bins = 21;
    int radius = 3;
    Grid grid(base, width, bins, radius);

    setting::grid::percent_water = 0;
    vector<Atom> a = {Atom({0, 0, 0}, 0, "C", "", 0)};
    grid.add(a);
    vector<Hetatom> water = grid.hydrate();

    IS_TRUE(water.size() == 6);
    if (water.size() == 6) { // avoid crashing if the above fails
        IS_TRUE(water[0].coords == Vector3({-6, 0, 0})); // (-2r, 0, 0)
        IS_TRUE(water[1].coords == Vector3({6, 0, 0})); // (2r, 0, 0)
        IS_TRUE(water[2].coords == Vector3({0, -6, 0})); // (0, -2r, 0)
        IS_TRUE(water[3].coords == Vector3({0, 6, 0})); // (0, 2r, 0)
        IS_TRUE(water[4].coords == Vector3({0, 0, -6})); // (0, 0, -2r)
        IS_TRUE(water[5].coords == Vector3({0, 0, 6})); // (0, 0, 2r)
    }
}

void test_width() {
    Vector3 base(-10, -10, -10);
    double width = 0.1;
    int bins = 210;
    int radius = 3;
    Grid grid(base, width, bins, radius);

    vector<Atom> a = {Atom({0, 0, 0}, 0, "C", "", 0)};
    grid.add(a);
    vector<vector<vector<char>>> &g = grid.grid;

    // check that it was placed correctly
    IS_TRUE(g[100][100][100] == 'A');

    // check water generation
    setting::grid::percent_water = 0;
    vector<Hetatom> water = grid.hydrate();
    IS_TRUE(water.size() == 6);
    if (water.size() == 6) { // avoid crashing if the above fails
        IS_TRUE(water[0].coords == Vector3({-6, 0, 0})); // (-2r, 0, 0)
        IS_TRUE(water[1].coords == Vector3({6, 0, 0})); // (2r, 0, 0)
        IS_TRUE(water[2].coords == Vector3({0, -6, 0})); // (0, -2r, 0)
        IS_TRUE(water[3].coords == Vector3({0, 6, 0})); // (0, 2r, 0)
        IS_TRUE(water[4].coords == Vector3({0, 0, -6})); // (0, 0, -2r)
        IS_TRUE(water[5].coords == Vector3({0, 0, 6})); // (0, 0, 2r)
    }

    // test bounds
    grid = Grid(base, width, bins, 3);

    a = {Atom({5, 0, -7}, 0, "C", "", 1), Atom({0, -5, 0}, 0, "C", "", 2)};
    grid.add(a);
    grid.expand_volume();
    g = grid.grid;

    vector<vector<int>> box = grid.bounding_box();
    IS_TRUE(box[0][0] == 100);
    IS_TRUE(box[0][1] == 150);
    IS_TRUE(box[1][0] == 50);
    IS_TRUE(box[1][1] == 100);
    IS_TRUE(box[2][0] == 30);
    IS_TRUE(box[2][1] == 100);
}

void test_remove() {
    Vector3 base(-10, -10, -10);
    double width = 1;
    int bins = 21;
    int radius = 3;
    Grid grid(base, width, bins, radius);

    Atom a1 = Atom({3, 0, 0}, 0, "C", "", 1);
    Atom a2 = Atom({0, 3, 0}, 0, "C", "", 2);
    vector<Atom> a = {a1, a2};
    grid.add(a);
    grid.expand_volume();
    grid.remove(a1);
    vector<vector<vector<char>>> &g = grid.grid;

    // check that the atom is removed
    IS_TRUE(g[13][10][10] == 0);
    IS_TRUE(g[13][10][9] == 0);
    IS_TRUE(g[13][11][10] == 0);
    IS_TRUE(g[13][11][9] == 0);

    // check the other is untouched
    IS_TRUE(g[10][13][10] == 'A');
    IS_TRUE(g[9][13][10] == 'a');
    IS_TRUE(g[10][12][10] == 'a');
    IS_TRUE(g[10][14][10] == 'a');
}

void test_find_free_locs(setting::grid::PlacementStrategyChoice ch) {
    Vector3 base(-10, -10, -10);
    double width = 1;
    vector<int> bins = {21, 21, 21};
    int radius = 3;
    Grid grid(base, width, bins, radius, radius, ch, setting::grid::csc);

    vector<Atom> a = {Atom({0, 0, 0}, 0, "C", "", 0)};
    grid.add(a);
    grid.expand_volume();

    vector<Hetatom> locs = grid.find_free_locs();
    IS_TRUE(locs.size() == 6);
    for (int i = 0; i < locs.size(); i++) {
        cout << format("%1%\t%2%\t%3%") % locs[i].coords.x % locs[i].coords.y % locs[i].coords.z << endl;
    }

    // since this needs to work with different placement strategies, we have to perform a more general check on the positions
    vector<Vector3> v = {{0, 0, 6}, {0, 0, -6}, {6, 0, 0}, {-6, 0, 0}, {0, 6, 0}, {0, -6, 0}};
    if (locs.size() >= 6) { // avoid crashing if the above fails
        for (const auto& l : locs) {
            bool found = false;
            for (const auto& p : v) {
                if (l.coords == p) {found = true;}
            }
            IS_TRUE(found);
        }
    }

}

int main(void)
{
    cout << "Summary of Grid testing:" << std::endl;
    test_grid_generation();
    cout << "1" << endl;
    test_simple_bounding_box();
    cout << "2" << endl;
    test_complex_bounding_box();
    cout << "3" << endl;
    test_find_free_locs(setting::grid::AxesStrategy);
    cout << "4" << endl;
    test_find_free_locs(setting::grid::RadialStrategy);
    cout << "5" << endl;
    test_hydrate();
    cout << "6" << endl;
    test_volume_expansion();
    cout << "7" << endl;
    test_width();
    cout << "8" << endl;
    test_remove();

    if (passed_all) {
        cout << "\033[1;32m" << "All Grid tests passed." << "\033[0m" << endl;
    } else {
        print_err("Some Grid tests failed.");
    }
    return 0;
}