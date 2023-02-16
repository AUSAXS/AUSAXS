#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <vector>
#include <string>
#include <fstream>

#include <data/Body.h>
#include <data/Protein.h>
#include <hydrate/Grid.h>
#include <hydrate/GridMember.h>
#include <utility/Settings.h>
#include <math/Vector3.h>
 
using std::vector;

// Debug class to expose the volume variable
class GridDebug : public Grid {
    using Grid::Grid;
    public: 
        auto unexpanded_volume() {return volume;}
        auto bounding_box_index() {return Grid::bounding_box_index();}
        auto to_bins(const Vector3<double> v) {return Grid::to_bins(v);}
        auto find_free_locs() {return Grid::find_free_locs();}
        static auto bounding_box(const vector<Atom>& atoms) {return Grid::bounding_box(atoms);}
};

TEST_CASE("grid_generation", "[grid]") {
    Axis3D axes(-10, 10, -10, 10, -10, 10, 20);
    int width = 1;
    Grid grid(axes, width);

    vector<Atom> a = {Atom({0, 0, 0}, 0, "C", "", 0)};
    grid.add(a);
    GridObj& g = grid.grid;

    // check that it was placed correctly in the grid
    REQUIRE(g.index(10, 10, 10) == GridObj::A_CENTER);
    REQUIRE(g.index(10, 10, 11) == GridObj::EMPTY);
    REQUIRE(g.index(10, 11, 11) == GridObj::EMPTY);
    REQUIRE(g.index(11, 10, 10) == GridObj::EMPTY);
    REQUIRE(g.index(9, 8, 14) == GridObj::EMPTY);
}

TEST_CASE("grid_member", "[grid]") {
    grid::GridMember<Atom> a(Atom({0.1, 0.2, 0.3}, 0, "C", "", 0), {1, 2, 3});
    auto b = a;

    CHECK(b == a);

    auto c(std::move(b));
    CHECK(c == a);
}

TEST_CASE("bounding_box", "[grid]") {
    Axis3D axes(-10, 10, -10, 10, -10, 10, 20);
    int width = 1;
    GridDebug grid(axes, width);

    SECTION("simple") {
        vector<Atom> a = {Atom({0, 0, 0}, 0, "C", "", 0)};
        grid.add(a);

        auto[min, max] = grid.bounding_box_index();
        CHECK(min[0] == 10);
        CHECK(max[0] == 11);
        CHECK(min[1] == 10);
        CHECK(max[1] == 11);
        CHECK(min[2] == 10);
        CHECK(max[2] == 11);
    }

    SECTION("complex") {
        vector<Atom> a = {Atom({5, 0, -7}, 0, "C", "", 1), Atom({0, -5, 0}, 0, "C", "", 2)};
        grid.add(a);
        grid.expand_volume();

        auto[min, max] = grid.bounding_box_index();
        CHECK(min[0] == 10);
        CHECK(max[0] == 16);
        CHECK(min[1] == 5);
        CHECK(max[1] == 11);
        CHECK(min[2] == 3);
        CHECK(max[2] == 11);
    }
}

TEST_CASE("volume_expansion", "[grid]") {
    Axis3D axes(-10, 10, -10, 10, -10, 10, 20);
    int width = 1;
    int radius = 3;
    Grid grid(axes, width, radius);

    vector<Atom> a = {Atom({0, 0, 0}, 0, "C", "", 0)};
    grid.add(a);
    GridObj &g = grid.grid;
    grid.expand_volume();

    REQUIRE(g.index(10, 10, 10) == GridObj::A_CENTER); // check that the center is still marked as capital 'A'

    // check that the x=13 plane looks like this: 
    // x x x
    // x o x
    // x x x
    CHECK(g.index(14, 10, 10) == GridObj::EMPTY); // check that it does not extend to the x=14 plane
    CHECK(g.index(13, 10, 10) == GridObj::A_AREA);

    CHECK(g.index(13, 9, 9) == GridObj::EMPTY);
    CHECK(g.index(13, 9, 10) == GridObj::EMPTY);
    CHECK(g.index(13, 9, 11) == GridObj::EMPTY);

    CHECK(g.index(13, 10, 9) == GridObj::EMPTY);
    CHECK(g.index(13, 11, 11) == GridObj::EMPTY);

    CHECK(g.index(13, 11, 9) == GridObj::EMPTY);
    CHECK(g.index(13, 11, 10) == GridObj::EMPTY);
    CHECK(g.index(13, 11, 11) == GridObj::EMPTY);

    // repeat with the x=7 plane
    CHECK(g.index(6, 10, 10) == GridObj::EMPTY); // check that it does not extend to the x=6 plane
    CHECK(g.index(7, 10, 10) == GridObj::A_AREA);

    CHECK(g.index(7, 9, 9) == GridObj::EMPTY);
    CHECK(g.index(7, 9, 10) == GridObj::EMPTY);
    CHECK(g.index(7, 9, 11) == GridObj::EMPTY);

    CHECK(g.index(7, 10, 9) == GridObj::EMPTY);
    CHECK(g.index(7, 11, 11) == GridObj::EMPTY);

    CHECK(g.index(7, 11, 9) == GridObj::EMPTY);
    CHECK(g.index(7, 11, 10) == GridObj::EMPTY);
    CHECK(g.index(7, 11, 11) == GridObj::EMPTY);

    // check some other points as well
    // x=10, z=10 line, from z=6 it looks like x o o o o o o o x
    CHECK(g.index(10, 6, 10) == GridObj::EMPTY);
    CHECK(g.index(10, 7, 10) == GridObj::A_AREA);
    CHECK(g.index(10, 8, 10) == GridObj::A_AREA);
    CHECK(g.index(10, 9, 10) == GridObj::A_AREA);
    CHECK(g.index(10, 10, 10) == GridObj::A_CENTER);
    CHECK(g.index(10, 11, 10) == GridObj::A_AREA);
    CHECK(g.index(10, 12, 10) == GridObj::A_AREA);
    CHECK(g.index(10, 13, 10) == GridObj::A_AREA);
    CHECK(g.index(10, 14, 10) == GridObj::EMPTY);

    // x=10, y=10 line, from y=6 it looks like x o o o o o o o x
    CHECK(g.index(10, 10, 6) == GridObj::EMPTY);
    CHECK(g.index(10, 10, 7) == GridObj::A_AREA);
    CHECK(g.index(10, 10, 8) == GridObj::A_AREA);
    CHECK(g.index(10, 10, 9) == GridObj::A_AREA);
    CHECK(g.index(10, 10, 10) == GridObj::A_CENTER);
    CHECK(g.index(10, 10, 11) == GridObj::A_AREA);
    CHECK(g.index(10, 10, 12) == GridObj::A_AREA);
    CHECK(g.index(10, 10, 13) == GridObj::A_AREA);
    CHECK(g.index(10, 10, 14) == GridObj::EMPTY);

    // some random points
    CHECK(g.index(9, 9, 9) == GridObj::A_AREA);
    CHECK(g.index(8, 8, 8) == GridObj::EMPTY);
    CHECK(g.index(13, 13, 13) == GridObj::EMPTY);
}

TEST_CASE("volume", "[grid]") {
    Axis3D axes(-10, 10, -10, 10, -10, 10, 20);
    int width = 1;
    int radius = 1;
    GridDebug grid(axes, width, radius);

    // cout << grid.get_volume() << endl;
    vector<Atom> a = {Atom({0, 0, 0}, 0, "C", "", 0)};
    vector<Water> w = {Water({2, 2, 2}, 0, "C", "", 0), Water({2, 2, 3}, 0, "C", "", 0)};
    grid.add(a);
    grid.add(w);
    REQUIRE(grid.unexpanded_volume() == 1); // atoms are added as point-particles, and only occupy one unit of space.

    grid.expand_volume();
    REQUIRE(grid.unexpanded_volume() == 7); // the radius is 1, so expanding the volume in a sphere results in one unit of volume added along each coordinate axis

    grid.add(Atom({0, 0, -1}, 0, "C", "", 0));
    grid.expand_volume();
    REQUIRE(grid.unexpanded_volume() == 12); // second atom is placed adjacent to the first one, so the volumes overlap. 
}

TEST_CASE("hydrate", "[grid],[files]") {
    setting::general::verbose = false;

    // check that all the expected hydration sites are found
    SECTION("correct placement") {
        Axis3D axes(-10, 10, -10, 10, -10, 10, 20);
        int width = 1;
        int radius = 3;
        Grid grid(axes, width, radius);

        // add a single atom to the grid, and hydrate it
        setting::grid::percent_water = 0;
        vector<Atom> a = {Atom({0, 0, 0}, 0, "C", "", 0)};
        grid.add(a);
        REQUIRE(grid.get_atoms().size() == 1); // check atoms

        vector<Water> water = grid.hydrate();

        REQUIRE(water.size() == 6); // check that the expected number of water molecules was placed
        CHECK(water[0].coords == Vector3{0, 0, 6});  // (0, 0,  2r)
        CHECK(water[1].coords == Vector3{0, 0, -6}); // (0, 0, -2r)
        CHECK(water[2].coords == Vector3{6, 0, 0});  // ( 2r, 0, 0)
        CHECK(water[3].coords == Vector3{-6, 0, 0}); // (-2r, 0, 0)
        CHECK(water[4].coords == Vector3{0, 6, 0});  // (0,  2r, 0)
        CHECK(water[5].coords == Vector3{0, -6, 0}); // (0, -2r, 0)

        REQUIRE(grid.get_atoms().size() == 1); // check that they are indeed registered as water molecules
    }

    // check that a hydration operation is reversible
    SECTION("reversible") {
        SECTION("LAR1-2") {
            // get the grid before hydrating it
            Protein protein("data/LAR1-2/LAR1-2.pdb");
            protein.clear_hydration();
            auto g1 = *protein.get_grid();

            // hydrate it and clear the hydration
            protein.generate_new_hydration();
            protein.clear_hydration();
            auto g2 = *protein.get_grid();

            // check sizes
            REQUIRE(g1.grid.xdim == g2.grid.xdim);
            REQUIRE(g1.grid.ydim == g2.grid.ydim);
            REQUIRE(g1.grid.zdim == g2.grid.zdim);

            // check that the grids are the same
            for (unsigned int i = 0; i < g1.grid.xdim; i++) {
                for (unsigned int j = 0; j < g1.grid.ydim; j++) {
                    for (unsigned int k = 0; k < g1.grid.zdim; k++) {
                        if (g1.grid.index(i, j, k) != g2.grid.index(i, j, k)) {
                            REQUIRE(g1.grid.index(i, j, k) != g2.grid.index(i, j, k));
                        }
                    }
                }
            }
            REQUIRE(g1.get_volume() == g2.get_volume());
        }
    }

    // check that a hydration operation produces consistent results
    SECTION("consistency") {
        Protein protein("data/LAR1-2/LAR1-2.pdb");
        protein.generate_new_hydration();
        auto h1 = protein.waters();
        auto a1 = protein.atoms();
        protein.generate_new_hydration();
        auto h2 = protein.waters();
        auto a2 = protein.atoms();

        // check that the hydration generation is deterministic
        for (unsigned int i = 0; i < h1.size(); i++) {
            REQUIRE(a1[i].coords == a2[i].coords);
        }
    }
}

TEST_CASE("width", "[grid]") {
    auto test_width_basics = [] (setting::grid::PlacementStrategy strategy) {
        setting::grid::placement_strategy = strategy;
        double width = 0.1;
        int radius = 3;
        Axis3D axes(-10, 10, -10, 10, -10, 10, width);
        Grid grid(axes, width, radius);

        vector<Atom> a = {Atom({0, 0, 0}, 0, "C", "", 0)};
        grid.add(a);
        GridObj& g = grid.grid;

        // check that it was placed correctly
        REQUIRE(g.index(100, 100, 100) == GridObj::A_CENTER);
        REQUIRE(grid.a_members.back().loc == Vector3(100, 100, 100));
        REQUIRE(grid.a_members.back().atom.coords == Vector3(0, 0, 0));

        // check water generation
        setting::grid::percent_water = 0;
        vector<Water> water = grid.hydrate();

        // since this needs to work with different placement strategies, we have to perform a more general check on the positions
        vector<Vector3<int>> v = {{0, 0, 6}, {0, 0, -6}, {6, 0, 0}, {-6, 0, 0}, {0, 6, 0}, {0, -6, 0}};
        REQUIRE(water.size() == 6);
        for (const auto& l : water) {
            bool found = false;
            for (const auto& p : v) {
                if (l.coords == p) {found = true;}
            }
            CHECK(found);
        }
    };

    auto test_width_bounds = [] (setting::grid::PlacementStrategy strategy) {
        setting::grid::placement_strategy = strategy;
        double width = 0.1;
        int radius = 3;
        Axis3D axes(-10, 10, -10, 10, -10, 10, width);
        GridDebug grid(axes, width, radius);

        vector<Atom> a = {Atom({5, 0, -7}, 0, "C", "", 1), Atom({0, -5, 0}, 0, "C", "", 2)};
        grid.add(a);
        grid.expand_volume();

        auto[min, max] = grid.bounding_box_index();
        REQUIRE(min[0] == 100);
        REQUIRE(max[0] == 151);
        REQUIRE(min[1] == 50);
        REQUIRE(max[1] == 101);
        REQUIRE(min[2] == 30);
        REQUIRE(max[2] == 101);
    };

    SECTION("RadialStrategy") {
        SECTION("basics") {
            test_width_basics(setting::grid::PlacementStrategy::RadialStrategy);
        }
        SECTION("bounds") {
            test_width_bounds(setting::grid::PlacementStrategy::RadialStrategy);
        }
    }
    SECTION("AxesStrategy") {
        SECTION("basics") {
            test_width_basics(setting::grid::PlacementStrategy::AxesStrategy);
        }
        SECTION("bounds") {
            test_width_bounds(setting::grid::PlacementStrategy::AxesStrategy);
        }
    }
    SECTION("JanStrategy") {
        SECTION("basics") {
            test_width_basics(setting::grid::PlacementStrategy::JanStrategy);
        }
        SECTION("bounds") {
            test_width_bounds(setting::grid::PlacementStrategy::JanStrategy);
        }
    }
}

TEST_CASE("add_remove", "[grid]") {
    Axis3D axes(-10, 10, -10, 10, -10, 10, 20);
    int width = 1;
    int radius = 3;
    GridDebug grid(axes, width, radius);

    // atoms
    Atom a1 = Atom({3, 0, 0}, 0, "C", "", 1);
    Atom a2 = Atom({0, 3, 0}, 0, "C", "", 2);
    Atom a3 = Atom({0, 0, 3}, 0, "C", "", 3);
    vector<Atom> a = {a1, a2, a3};

    // waters
    Water w1 = Water::create_new_water(Vector3<double>({0, 0, -3}));
    Water w2 = Water::create_new_water(Vector3<double>({0, -3, 0}));
    Water w3 = Water::create_new_water(Vector3<double>({-3, 0, 0}));
    vector<Water> w = {w1, w2, w3};

    SECTION("add") {
        // add atoms
        grid.add(a);    
        vector<Atom> ga = grid.get_atoms();
        REQUIRE(grid.a_members.size() == 3); // check actual data
        REQUIRE(ga.size() >= 3);
        for (int i = 0; i < 3; i++) {
            REQUIRE(ga[i] == a[i]); // check equality with what we added
        }

        // add waters
        grid.add(w);
        vector<Water> wa = grid.get_waters();
        REQUIRE(grid.a_members.size() == 3); // check actual data
        REQUIRE(grid.w_members.size() == 3); // check actual data
        REQUIRE(wa.size() >= 3);
        for (int i = 0; i < 3; i++) {
            REQUIRE(wa[i] == w[i]); // check equality with what we added
        }
    }

    SECTION("remove") {
        grid.add(a);
        grid.add(w);

        grid.remove(a2);
        grid.remove(w3);
        grid.remove(w1);
        vector<Atom> ga = grid.get_atoms();
        vector<Water> wa = grid.get_waters();

        // check sizes
        REQUIRE(ga.size() == 2);
        REQUIRE(wa.size() == 1);

        // check remaining atoms
        REQUIRE(ga[0] == a1);
        REQUIRE(ga[1] == a3);
        REQUIRE(wa[0] == w2);

        // check old centers
        GridObj &g = grid.grid;
        auto loc_a2 = grid.to_bins(a2.coords);
        auto loc_w1 = grid.to_bins(w1.coords);
        auto loc_w3 = grid.to_bins(w3.coords);
        REQUIRE(g.index(loc_a2) == GridObj::EMPTY);
        REQUIRE(g.index(loc_w1) == GridObj::EMPTY);
        REQUIRE(g.index(loc_w3) == GridObj::EMPTY);
    }

    SECTION("remove_vector") {
        grid.add(a);
        grid.add(w);

        // remove a list of hetatoms
        vector<Water> remove_water = {w1, w3};
        grid.remove(remove_water);

        // check the remaining one
        vector<Atom> ga = grid.get_atoms();
        vector<Water> wa = grid.get_waters();
        REQUIRE(wa.size() == 1);
        REQUIRE(ga.size() == 3);
        REQUIRE(wa[0] == w2);
    }

    SECTION("remove index") {
        Atom a4 = Atom({1, 0, 0}, 0, "C", "", 4);
        Atom a5 = Atom({0, 1, 0}, 0, "C", "", 5);
        Atom a6 = Atom({0, 0, 1}, 0, "C", "", 6);
        Atom a7 = Atom({2, 0, 0}, 0, "C", "", 7);
        Atom a8 = Atom({0, 2, 0}, 0, "C", "", 8);
        Atom a9 = Atom({0, 0, 2}, 0, "C", "", 9);
        a = {a1, a2, a3, a4, a5, a6, a7, a8, a9}; 

        grid.add(a);
        std::vector<bool> remove = {false, true, false, false, false, false, true, true, true};
        grid.remove(remove);

        vector<Atom> ga = grid.get_atoms();
        REQUIRE(ga.size() == 5);
        REQUIRE(ga[0] == a1);
        REQUIRE(ga[1] == a3);
        REQUIRE(ga[2] == a4);
        REQUIRE(ga[3] == a5);
        REQUIRE(ga[4] == a6);        
    }

    SECTION("clear_waters") {
        grid.add(a);
        grid.add(w);

        grid.clear_waters();
        REQUIRE(grid.a_members.size() == 3);
        REQUIRE(grid.w_members.size() == 0);
    }
}

TEST_CASE("correct_volume", "[grid]") {
    Axis3D axes(-10, 10, -10, 10, -10, 10, 20);
    int width = 1;
    int radius = 10; // we use a very large radius so the atoms will heavily overlap
    GridDebug grid(axes, width, radius);

    // atoms
    Atom a1 = Atom({3, 0, 0}, 0, "C", "", 1);
    Atom a2 = Atom({0, 3, 0}, 0, "C", "", 2);
    Atom a3 = Atom({0, 0, 3}, 0, "C", "", 3);
    vector<Atom> a = {a1, a2, a3};

    // waters
    Water w1 = Water::create_new_water(Vector3<double>({0, 0, -3}));
    Water w2 = Water::create_new_water(Vector3<double>({0, -3, 0}));
    Water w3 = Water::create_new_water(Vector3<double>({-3, 0, 0}));
    vector<Water> w = {w1, w2, w3};

    // single non-overlapping
    REQUIRE(grid.unexpanded_volume() == 0);
    grid.add(a1);
    REQUIRE(grid.unexpanded_volume() != 0);
    grid.remove(a1);
    REQUIRE(grid.unexpanded_volume() == 0);

    grid.add(w1);
    REQUIRE(grid.unexpanded_volume() == 0);
    grid.remove(w1);
    REQUIRE(grid.unexpanded_volume() == 0);

    // multiple overlapping
    grid.add(a);
    grid.add(w);
    REQUIRE(grid.unexpanded_volume() != 0);
    grid.remove(a);
    grid.remove(w);
    REQUIRE(grid.unexpanded_volume() == 0);
}

TEST_CASE("find_free_locs", "[grid]") {
    auto test_func = [] (setting::grid::PlacementStrategy ch) {
        Axis3D axes(-10, 10, -10, 10, -10, 10, 20);
        int width = 1;
        int radius = 3;
        GridDebug grid(axes, width, radius, radius, ch, setting::grid::culling_strategy);

        vector<Atom> a = {Atom({0, 0, 0}, 0, "C", "", 0)};
        grid.add(a);
        grid.expand_volume();

        vector<grid::GridMember<Water>> locs = grid.find_free_locs();
        REQUIRE(locs.size() == 6);

        // since this needs to work with different placement strategies, we have to perform a more general check on the positions
        vector<Vector3<int>> v = {{0, 0, 6}, {0, 0, -6}, {6, 0, 0}, {-6, 0, 0}, {0, 6, 0}, {0, -6, 0}};
        for (const auto& l : locs) {
            bool found = false;
            for (const auto& p : v) {
                if (l.atom.coords == p) {found = true;}
            }
            REQUIRE(found);
        }
    };

    SECTION("RadialStrategy") {
        test_func(setting::grid::PlacementStrategy::RadialStrategy);
    }
    SECTION("JanStrategy") {
        test_func(setting::grid::PlacementStrategy::JanStrategy);
    }
    SECTION("AxesStrategy") {
        test_func(setting::grid::PlacementStrategy::AxesStrategy);
    }
}

// Test that expansion and deflation completely cancels each other. 
TEST_CASE("volume_deflation", "[grid]") {
    Axis3D axes(-10, 10, -10, 10, -10, 10, 20);
    int width = 1;
    int radius = 3;
    GridDebug grid(axes, width, radius);

    vector<Atom> a = {Atom({3, 0, 0}, 0, "C", "", 1), Atom({0, 3, 0}, 0, "C", "", 2)};
    grid.add(a);
    REQUIRE(grid.unexpanded_volume() == 2);

    grid.expand_volume();
    grid.deflate_volume();
    GridObj &g = grid.grid;
    REQUIRE(grid.unexpanded_volume() == 2);

    auto bins = grid.get_bins();
    for (int i = 0; i < bins.x(); i++) {
        for (int j = 0; j < bins.y(); j++) {
            for (int k = 0; k < bins.z(); k++) {
                if (i == 10 && j == 13 && k == 10, false) {continue;} // center of the first atom
                if (i == 13 && j == 10 && k == 10, false) {continue;} // center of the second atom
                if (g.index(i, j, k) != GridObj::EMPTY) {
                    REQUIRE(false);
                }
            }
        }
    }
}

TEST_CASE("grid_cubic", "[grid]") {
    setting::grid::cubic = true;

    SECTION("largest x") {
        Axis3D axes(-10, 10, -1, 1, -1, 1, 20);
        int width = 1;
        int radius = 3;

        Grid grid(axes, width, radius);
        auto gaxes = grid.get_axes();
        CHECK(gaxes.x == axes.x);
        CHECK(gaxes.x == gaxes.y);
        CHECK(gaxes.x == gaxes.z);
    }

    SECTION("largest y") {
        Axis3D axes(-1, 1, -10, 10, -1, 1, 20);
        int width = 1;
        int radius = 3;

        Grid grid(axes, width, radius);
        auto gaxes = grid.get_axes();
        CHECK(gaxes.y == axes.y);
        CHECK(gaxes.y == gaxes.y);
        CHECK(gaxes.y == gaxes.z);
    }

    SECTION("largest x") {
        Axis3D axes(-1, 1, -1, 1, -10, 10, 20);
        int width = 1;
        int radius = 3;

        Grid grid(axes, width, radius);
        auto gaxes = grid.get_axes();
        CHECK(gaxes.z == axes.z);
        CHECK(gaxes.z == gaxes.y);
        CHECK(gaxes.z == gaxes.z);
    }

    setting::grid::cubic = false;
}

TEST_CASE("grid_space_saving_constructor", "[grid]") {
    vector<Atom> atoms = {Atom({5, 0, -7}, 0, "C", "", 1), Atom({0, -5, 0}, 0, "C", "", 2), Atom({1, 1, 1}, 0, "C", "", 2)};

    SECTION("bounding box") {
        // check that bounding_box works
        auto[min, max] = GridDebug::bounding_box(atoms);
        CHECK(min.x() == 0);
        CHECK(min.y() == -5);
        CHECK(min.z() == -7);
        CHECK(max.x() == 5);
        CHECK(max.y() == 1);
        CHECK(max.z() == 1);
    }

    auto func = [] (const Grid& grid) {
        Axis3D axes = grid.get_axes();
        CHECK(axes.x.min ==  0 - 5*setting::grid::scaling*0.5 - setting::grid::width);
        CHECK(axes.y.min == -5 - 6*setting::grid::scaling*0.5 - setting::grid::width);
        CHECK(axes.z.min == -7 - 8*setting::grid::scaling*0.5 - setting::grid::width);
        CHECK(axes.x.max ==  5 + 5*setting::grid::scaling*0.5 + setting::grid::width);
        CHECK(axes.y.max ==  1 + 6*setting::grid::scaling*0.5 + setting::grid::width);
        CHECK(axes.z.max ==  1 + 8*setting::grid::scaling*0.5 + setting::grid::width);

        // check that we're not using a ton of unnecessary bins
        CHECK(axes.x.bins < 20);
        CHECK(axes.y.bins < 20);
        CHECK(axes.z.bins < 20);

        // check that this is reflected in the grid itself
        const GridObj& g = grid.grid;
        CHECK(g.xdim < 20);
        CHECK(g.ydim < 20);
        CHECK(g.zdim < 20);
    };

    SECTION("single body") {
        // check that the grid constructor works as expected
        Grid grid(atoms);
        func(grid);
    }

    SECTION("multiple bodies") {
        vector<Body> bodies; 
        for (const auto& a : atoms) {bodies.push_back(Body({a}));}

        Grid grid(bodies);
        func(grid);
    }
}

TEST_CASE("grid_copy", "[grid]") {
    Axis3D axes(-10, 10, -10, 10, -10, 10, 20);
    Grid grid1(axes, 1);
    grid1.add(Atom({0, 0, 0}, 0, "C", "", 0));
    grid1.hydrate();

    // copy
    Grid grid2 = grid1.copy();
    REQUIRE(grid2 == grid1);

    // assignment
    Grid grid3(axes, 1);
    grid3 = grid1;
    REQUIRE(grid3 == grid1);
}