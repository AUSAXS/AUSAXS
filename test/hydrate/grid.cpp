#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <vector>
#include <string>

#include <data/Body.h>
#include <data/Molecule.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <grid/Grid.h>
#include <grid/detail/GridMember.h>
#include <hydrate/generation/HydrationStrategy.h>
#include <hydrate/generation/HydrationFactory.h>
#include <hydrate/culling/CullingStrategy.h>
#include <settings/All.h>
#include <math/Vector3.h>
#include <constants/Constants.h>
 
using std::vector;
using namespace grid;
using namespace grid::detail;
using namespace data;
using namespace data::record;

// Debug class to expose the volume variable
class GridDebug : public grid::Grid {
    public: 
        GridDebug(Limit3D axes) : Grid(axes) {}

		double get_atomic_radius(constants::atom_t) const override {return ra;}
		double get_hydration_radius() const override {return rh;}
        void set_atomic_radius(double ra) {this->ra = ra;}
        void set_hydration_radius(double rh) {this->rh = rh;}
        auto get_volume_without_expanding() {return volume;}
        auto bounding_box_index() {return Grid::bounding_box_index();}
        auto to_bins(const Vector3<double> v) {return Grid::to_bins(v);}
        auto get_member_atoms() {return a_members;}
        auto get_member_waters() {return w_members;}
        static auto bounding_box(const vector<Atom>& atoms) {return Grid::bounding_box(atoms);}
    
    private:
        double ra = 0, rh = 0;
};

TEST_CASE("Grid::Grid") {
    settings::grid::width = 1;

    SECTION("Limit3D&") {
        Limit3D axes(-10, 10, -10, 10, -10, 10);
        Grid grid(axes);
        CHECK(grid.a_members.empty());
        CHECK(grid.w_members.empty());
        CHECK(grid.get_atoms().empty());
        CHECK(grid.get_volume() == 0);
        CHECK(grid.get_width() == settings::grid::width);
        CHECK(grid.get_axes() == Axis3D(axes, settings::grid::width));
    }

    SECTION("vector<Atom>&") {
        std::vector<Atom> atoms = {Atom({0, 0, 0}, 0,  constants::atom_t::C, "", 0), Atom({1, 1, 1}, 0,  constants::atom_t::C, "", 0), Atom({2, 2, 2}, 0,  constants::atom_t::C, "", 0)};
        Grid grid(atoms);
        CHECK(grid.a_members.size() == 3);
        CHECK(grid.w_members.empty());
        CHECK(grid.get_atoms() == atoms);
        CHECK(grid.get_volume() != 0);
        CHECK(grid.get_width() == settings::grid::width);
    }

    SECTION("std::vector<Body>&") {
        std::vector<Body> bodies = {Body({Atom({0, 0, 0}, 0,  constants::atom_t::C, "", 0), Atom({1, 1, 1}, 0,  constants::atom_t::C, "", 0), Atom({2, 2, 2}, 0,  constants::atom_t::C, "", 0)})};
        Grid grid(bodies);
        CHECK(grid.a_members.size() == 3);
        CHECK(grid.w_members.empty());
        CHECK(grid.get_atoms() == bodies[0].get_atoms());
        CHECK(grid.get_volume() != 0);
        CHECK(grid.get_width() == settings::grid::width);
    }

    SECTION("space_saving_constructor") {
        vector<Atom> atoms = {Atom({5, 0, -7}, 0,  constants::atom_t::C, "", 1), Atom({0, -5, 0}, 0,  constants::atom_t::C, "", 2), Atom({1, 1, 1}, 0,  constants::atom_t::C, "", 2)};

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
            CHECK(axes.x.min == std::floor( 0 - 5*settings::grid::scaling*0.5 - settings::grid::width));
            CHECK(axes.y.min == std::floor(-5 - 6*settings::grid::scaling*0.5 - settings::grid::width));
            CHECK(axes.z.min == std::floor(-7 - 8*settings::grid::scaling*0.5 - settings::grid::width));
            CHECK(axes.x.max ==  std::ceil(5 + 5*settings::grid::scaling*0.5 + settings::grid::width));
            CHECK(axes.y.max ==  std::ceil(1 + 6*settings::grid::scaling*0.5 + settings::grid::width));
            CHECK(axes.z.max ==  std::ceil(1 + 8*settings::grid::scaling*0.5 + settings::grid::width));

            // check that we're not using a ton of unnecessary bins
            CHECK(axes.x.bins < 20);
            CHECK(axes.y.bins < 20);
            CHECK(axes.z.bins < 20);

            // check that this is reflected in the grid itself
            const GridObj& g = grid.grid;
            CHECK(g.size_x() < 20);
            CHECK(g.size_y() < 20);
            CHECK(g.size_z() < 20);
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
}

TEST_CASE("generation") {
    Limit3D axes(-10, 10, -10, 10, -10, 10);
    settings::grid::width = 1;
    Grid grid(axes);

    vector<Atom> a = {Atom({0, 0, 0}, 0,  constants::atom_t::C, "", 0)};
    grid.add(a);
    GridObj& g = grid.grid;

    // check that it was placed correctly in the grid
    REQUIRE(g.index(10, 10, 10) == A_CENTER);
    REQUIRE(g.index(10, 10, 11) == EMPTY);
    REQUIRE(g.index(10, 11, 11) == EMPTY);
    REQUIRE(g.index(11, 10, 10) == EMPTY);
    REQUIRE(g.index(9, 8, 14) == EMPTY);
}

TEST_CASE("GridMember") {
    grid::GridMember<Atom> a(Atom({0.1, 0.2, 0.3}, 0,  constants::atom_t::C, "", 0), {1, 2, 3});

    auto b = a;
    CHECK(b == a);

    auto c(std::move(b));
    CHECK(c == a);
}

TEST_CASE("Grid::bounding_box") {
    Limit3D axes(-10, 10, -10, 10, -10, 10);
    settings::grid::width = 1;
    GridDebug grid(axes);

    SECTION("simple") {
        vector<Atom> a = {Atom({0, 0, 0}, 0,  constants::atom_t::C, "", 0)};
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
        vector<Atom> a = {Atom({5, 0, -7}, 0,  constants::atom_t::C, "", 1), Atom({0, -5, 0}, 0,  constants::atom_t::C, "", 2)};
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

TEST_CASE("Grid::get_volume") {
    Limit3D lims(-10, 10, -10, 10, -10, 10);
    settings::grid::width = 1e-1;
    Grid grid(lims);

    std::vector<Atom> a = {Atom({0, 0, 0}, 0,  constants::atom_t::C, "", 0)};
    grid.add(a);
    grid.expand_volume();
    GridObj &g = grid.grid;

    unsigned int count = 0;
    auto axes = grid.get_axes();
    for (unsigned int i = 0; i < axes.x.bins; i++) {
        for (unsigned int j = 0; j < axes.y.bins; j++) {
            for (unsigned int k = 0; k < axes.z.bins; k++) {
                if (g.index(i, j, k) != EMPTY) {
                    count++;
                }
            }
        }
    }

    REQUIRE(count != 0);
    CHECK(grid.get_volume() == count*std::pow(settings::grid::width, 3));
}

TEST_CASE("Grid::width") {
    // check that the grid width is actually used
    settings::grid::width = GENERATE(0.25, 0.5, 1, 2);

    Limit3D lims(-10, 10, -10, 10, -10, 10);
    Grid grid(lims);

    auto axes = grid.get_axes();
    for (unsigned int i = 0; i < axes.x.bins-1; ++i) {
        CHECK(grid.to_xyz(i, 0, 0).distance(grid.to_xyz(i+1, 0, 0)) == settings::grid::width);
        CHECK(grid.to_xyz(0, i, 0).distance(grid.to_xyz(0, i+1, 0)) == settings::grid::width);
        CHECK(grid.to_xyz(0, 0, i).distance(grid.to_xyz(0, 0, i+1)) == settings::grid::width);
    }

    auto center = grid.get_center();
    CHECK(grid.to_xyz(center) == Vector3<double>{0, 0, 0});
    CHECK(grid.to_xyz(center+Vector3<double>{1, 0, 0}) == Vector3<double>{settings::grid::width, 0, 0});
    CHECK(grid.to_xyz(center+Vector3<double>{0, 1, 0}) == Vector3<double>{0, settings::grid::width, 0});
    CHECK(grid.to_xyz(center+Vector3<double>{0, 0, 1}) == Vector3<double>{0, 0, settings::grid::width});
}

TEST_CASE("Grid::expand_volume") {
    settings::molecule::use_effective_charge = false;
    Limit3D lims(-10, 10, -10, 10, -10, 10);
    settings::grid::width = 1;
    Grid grid(lims);
    constants::radius::set_dummy_radius(3+1e-6);

    vector<Atom> a = {Atom({0, 0, 0}, 0,  constants::atom_t::dummy, "", 0)};
    grid.add(a);
    grid.expand_volume();

    auto axes = grid.get_axes();
    for (unsigned int i = 0; i < axes.x.bins; ++i) {
        for (unsigned int j = 0; j < axes.y.bins; ++j) {
            for (unsigned int k = 0; k < axes.z.bins; ++k) {
                auto& bin = grid.grid.index(i, j, k);
                if (grid.to_xyz(i, j, k).norm() <= 3) {
                    if (i == 10 && j == 10 && k == 10) {continue;}
                    CHECK(bin == A_AREA);
                } else {
                    CHECK(bin == EMPTY);
                }
            }
        }
    }
    CHECK(grid.grid.index(10, 10, 10) == A_CENTER);
}

TEST_CASE("volume") {
    Limit3D axes(-10, 10, -10, 10, -10, 10);
    settings::grid::width = 1;
    GridDebug grid(axes);
    grid.set_atomic_radius(1);

    // cout << grid.get_volume() << endl;
    vector<Atom> a = {Atom({0, 0, 0}, 0,  constants::atom_t::dummy, "", 0)};
    vector<Water> w = {Water({2, 2, 2}, 0,  constants::atom_t::dummy, "", 0), Water({2, 2, 3}, 0,  constants::atom_t::dummy, "", 0)};
    grid.add(a);
    grid.add(w);
    REQUIRE(grid.get_volume_without_expanding() == 1); // atoms are added as point-particles, and only occupy one unit of space.

    SECTION("") {
        settings::grid::rvol = 1;
        REQUIRE(grid.get_volume() == 7); // the radius is 1, so expanding the volume in a sphere results in one unit of volume added along each coordinate axis

        grid.add(Atom({0, 0, -1}, 0,  constants::atom_t::C, "", 0));
        REQUIRE(grid.get_volume() == 12); // second atom is placed adjacent to the first one, so the volumes overlap. 
    }

    SECTION("") {
        settings::grid::rvol = 2;
        grid.expand_volume();
        REQUIRE(grid.get_volume() == 33);
    }
}

TEST_CASE("Grid::hydrate") {
    settings::general::verbose = false;
    settings::molecule::center = false;
    settings::hydrate::hydration_strategy = GENERATE(settings::hydrate::HydrationStrategy::AxesStrategy, settings::hydrate::HydrationStrategy::RadialStrategy, settings::hydrate::HydrationStrategy::JanStrategy);

    // check that all the expected hydration sites are found
    SECTION("correct placement " + std::to_string(static_cast<int>(settings::hydrate::hydration_strategy))) {
        Limit3D axes(-10, 10, -10, 10, -10, 10);
        settings::grid::width = 1;
        auto grid = std::make_unique<GridDebug>(axes);
        grid->set_atomic_radius(3);
        grid->set_hydration_radius(3);

        // add a single atom to the grid, and hydrate it
        settings::grid::water_scaling = 0;
        vector<Atom> a = {Atom({0, 0, 0}, 0,  constants::atom_t::dummy, "", 0)};
        grid->add(a);
        REQUIRE(grid->get_atoms().size() == 1); // check atoms

        data::Molecule protein(a);
        protein.set_grid(std::move(grid));
        protein.generate_new_hydration();
        auto water = protein.get_waters();

        REQUIRE(6 <= water.size()); // check that the expected number of water molecules was placed
        int count = 0;
        for (const auto& w : water) {
            if (w.coords == Vector3{ 0,  0,  6}) {++count; continue;} // ( 0,   0,  2r)
            if (w.coords == Vector3{ 0,  0, -6}) {++count; continue;} // ( 0,   0, -2r)
            if (w.coords == Vector3{ 6,  0,  0}) {++count; continue;} // ( 2r,  0,   0)
            if (w.coords == Vector3{-6,  0,  0}) {++count; continue;} // (-2r,  0,   0)
            if (w.coords == Vector3{ 0,  6,  0}) {++count; continue;} // ( 0,  2r,   0)
            if (w.coords == Vector3{ 0, -6,  0}) {++count; continue;} // ( 0, -2r,   0)
            std::cout << w.coords << std::endl;
        }
        REQUIRE(count == 6);
        REQUIRE(protein.get_grid()->get_atoms().size() == 1); // check that they are indeed registered as water molecules
    }

    // check that a hydration operation is reversible
    SECTION("reversible") {
        SECTION("LAR1-2") {
            // get the grid before hydrating it
            Molecule protein("test/files/LAR1-2.pdb");
            protein.clear_hydration();
            auto g1 = *protein.get_grid();

            // hydrate it and clear the hydration
            protein.generate_new_hydration();
            protein.clear_hydration();
            auto g2 = *protein.get_grid();

            // check sizes
            REQUIRE(g1.grid.size_x() == g2.grid.size_x());
            REQUIRE(g1.grid.size_y() == g2.grid.size_y());
            REQUIRE(g1.grid.size_z() == g2.grid.size_z());

            // check that the grids are the same
            for (unsigned int i = 0; i < g1.grid.size_x(); i++) {
                for (unsigned int j = 0; j < g1.grid.size_y(); j++) {
                    for (unsigned int k = 0; k < g1.grid.size_z(); k++) {
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
        settings::hydrate::culling_strategy = settings::hydrate::CullingStrategy::CounterStrategy;
        settings::hydrate::hydration_strategy = GENERATE(settings::hydrate::HydrationStrategy::AxesStrategy, settings::hydrate::HydrationStrategy::RadialStrategy, settings::hydrate::HydrationStrategy::JanStrategy);
        Molecule protein("test/files/LAR1-2.pdb");
        protein.get_grid()->force_expand_volume();

        protein.generate_new_hydration();
        auto h1 = protein.get_waters();

        protein.generate_new_hydration();
        auto h2 = protein.get_waters();

        // check that the hydration generation is deterministic
        REQUIRE(h1.size() == h2.size());
        for (unsigned int i = 0; i < h1.size(); ++i) {
            REQUIRE(h1[i].coords == h2[i].coords);
        }
    }
}

TEST_CASE("Grid: using different widths") {
    settings::general::verbose = false;
    auto test_width_basics = [] (settings::hydrate::HydrationStrategy strategy) {
        settings::hydrate::hydration_strategy = strategy;
        settings::grid::width = 0.1;
        vector<Atom> a = {Atom({0, 0, 0}, 0,  constants::atom_t::C, "LYS", 0)};
        data::Molecule protein(a);
        {
            Limit3D axes(-10, 10, -10, 10, -10, 10);
            auto grid = std::make_unique<GridDebug>(axes);
            grid->set_atomic_radius(3);
            grid->set_hydration_radius(3);

            grid->add(a);
            GridObj& g = grid->grid;

            // check that it was placed correctly
            REQUIRE(g.index(100, 100, 100) == A_CENTER);
            REQUIRE(grid->a_members.back().get_bin_loc() == Vector3(100, 100, 100));
            REQUIRE(grid->a_members.back().get_atom().get_coordinates() == Vector3(0, 0, 0));

            // check water generation
            settings::grid::water_scaling = 0;

            protein.set_grid(std::move(grid));
            protein.generate_new_hydration();
        }
        vector<Water> water = protein.get_waters();

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

    auto test_width_bounds = [] (settings::hydrate::HydrationStrategy strategy) {
        settings::hydrate::hydration_strategy = strategy;
        settings::grid::width = 0.1;
        Limit3D axes(-10, 10, -10, 10, -10, 10);
        GridDebug grid(axes);
        grid.set_atomic_radius(3);
        grid.set_hydration_radius(3);

        vector<Atom> a = {Atom({5, 0, -7}, 0,  constants::atom_t::C, "", 1), Atom({0, -5, 0}, 0,  constants::atom_t::C, "", 2)};
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
            test_width_basics(settings::hydrate::HydrationStrategy::RadialStrategy);
        }
        SECTION("bounds") {
            test_width_bounds(settings::hydrate::HydrationStrategy::RadialStrategy);
        }
    }
    SECTION("AxesStrategy") {
        SECTION("basics") {
            test_width_basics(settings::hydrate::HydrationStrategy::AxesStrategy);
        }
        SECTION("bounds") {
            test_width_bounds(settings::hydrate::HydrationStrategy::AxesStrategy);
        }
    }
    SECTION("JanStrategy") {
        SECTION("basics") {
            test_width_basics(settings::hydrate::HydrationStrategy::JanStrategy);
        }
        SECTION("bounds") {
            test_width_bounds(settings::hydrate::HydrationStrategy::JanStrategy);
        }
    }
}

TEST_CASE("Grid: add and remove") {
    Limit3D axes(-10, 10, -10, 10, -10, 10);
    settings::grid::width = 1;
    GridDebug grid(axes);
    grid.set_atomic_radius(3);
    grid.set_hydration_radius(3);

    // atoms
    Atom a1 = Atom({3, 0, 0}, 0,  constants::atom_t::C, "", 1);
    Atom a2 = Atom({0, 3, 0}, 0,  constants::atom_t::C, "", 2);
    Atom a3 = Atom({0, 0, 3}, 0,  constants::atom_t::C, "", 3);
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
        REQUIRE(g.index(loc_a2) == EMPTY);
        REQUIRE(g.index(loc_w1) == EMPTY);
        REQUIRE(g.index(loc_w3) == EMPTY);
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
        Atom a4 = Atom({1, 0, 0}, 0,  constants::atom_t::C, "", 4);
        Atom a5 = Atom({0, 1, 0}, 0,  constants::atom_t::C, "", 5);
        Atom a6 = Atom({0, 0, 1}, 0,  constants::atom_t::C, "", 6);
        Atom a7 = Atom({2, 0, 0}, 0,  constants::atom_t::C, "", 7);
        Atom a8 = Atom({0, 2, 0}, 0,  constants::atom_t::C, "", 8);
        Atom a9 = Atom({0, 0, 2}, 0,  constants::atom_t::C, "", 9);
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

TEST_CASE("Grid: correct_volume") {
    Limit3D axes(-10, 10, -10, 10, -10, 10);
    settings::grid::width = 1;
    GridDebug grid(axes);
    grid.set_atomic_radius(10);     // heavy overlap
    grid.set_hydration_radius(10);

    // atoms
    Atom a1 = Atom({3, 0, 0}, 0,  constants::atom_t::C, "", 1);
    Atom a2 = Atom({0, 3, 0}, 0,  constants::atom_t::C, "", 2);
    Atom a3 = Atom({0, 0, 3}, 0,  constants::atom_t::C, "", 3);
    vector<Atom> a = {a1, a2, a3};

    // waters
    Water w1 = Water::create_new_water(Vector3<double>({0, 0, -3}));
    Water w2 = Water::create_new_water(Vector3<double>({0, -3, 0}));
    Water w3 = Water::create_new_water(Vector3<double>({-3, 0, 0}));
    vector<Water> w = {w1, w2, w3};

    // single non-overlapping
    REQUIRE(grid.get_volume_without_expanding() == 0);
    grid.add(a1);
    REQUIRE(grid.get_volume_without_expanding() != 0);
    grid.remove(a1);
    REQUIRE(grid.get_volume_without_expanding() == 0);

    grid.add(w1);
    REQUIRE(grid.get_volume_without_expanding() == 0);
    grid.remove(w1);
    REQUIRE(grid.get_volume_without_expanding() == 0);

    // multiple overlapping
    grid.add(a);
    grid.add(w);
    REQUIRE(grid.get_volume_without_expanding() != 0);
    grid.remove(a);
    grid.remove(w);
    REQUIRE(grid.get_volume_without_expanding() == 0);
}

TEST_CASE("Grid::find_free_locs") {
    settings::general::verbose = false;
    auto test_func = [] (settings::hydrate::HydrationStrategy ch) {
        settings::hydrate::hydration_strategy = ch;
        settings::grid::width = 1;

        Limit3D axes(-10, 10, -10, 10, -10, 10);
        auto grid = std::make_unique<GridDebug>(axes);
        grid->set_atomic_radius(3);
        grid->set_hydration_radius(3);

        settings::grid::water_scaling = 0;
        vector<Atom> a = {Atom({0, 0, 0}, 0,  constants::atom_t::C, "LYS", 0)};
        grid->add(a);
        grid->expand_volume();

        data::Molecule protein(a);
        protein.set_grid(std::move(grid));
        protein.generate_new_hydration();

        std::list<grid::GridMember<Water>> locs = static_cast<GridDebug*>(protein.get_grid())->get_member_waters();
        REQUIRE(locs.size() == 6);

        // since this needs to work with different placement strategies, we have to perform a more general check on the positions
        vector<Vector3<int>> v = {{0, 0, 6}, {0, 0, -6}, {6, 0, 0}, {-6, 0, 0}, {0, 6, 0}, {0, -6, 0}};
        for (const auto& l : locs) {
            bool found = false;
            for (const auto& p : v) {
                if (l.get_atom().get_coordinates() == p) {found = true;}
            }
            REQUIRE(found);
        }
    };

    SECTION("RadialStrategy") {
        test_func(settings::hydrate::HydrationStrategy::RadialStrategy);
    }
    SECTION("JanStrategy") {
        test_func(settings::hydrate::HydrationStrategy::JanStrategy);
    }
    SECTION("AxesStrategy") {
        test_func(settings::hydrate::HydrationStrategy::AxesStrategy);
    }
}

// Test that expansion and deflation completely cancels each other. 
TEST_CASE("Grid::deflate_volume") {
    settings::general::verbose = false;

    Limit3D axes(-10, 10, -10, 10, -10, 10);
    settings::grid::width = 1;
    auto grid = std::make_unique<GridDebug>(axes);
    grid->set_atomic_radius(3);
    grid->set_hydration_radius(3);

    vector<Atom> a = {Atom({3, 0, 0}, 0,  constants::atom_t::C, "", 1), Atom({0, 3, 0}, 0,  constants::atom_t::C, "", 2)};
    grid->add(a);
    REQUIRE(grid->get_volume_without_expanding() == 2);

    grid->expand_volume();
    grid->deflate_volume();
    GridObj &g = grid->grid;
    REQUIRE(grid->get_volume_without_expanding() == 2);

    auto bins = grid->get_bins();
    for (int i = 0; i < bins.x(); i++) {
        for (int j = 0; j < bins.y(); j++) {
            for (int k = 0; k < bins.z(); k++) {
                if (i == 10 && j == 13 && k == 10) {continue;} // center of the first atom
                if (i == 13 && j == 10 && k == 10) {continue;} // center of the second atom
                if (g.index(i, j, k) != EMPTY) {
                    REQUIRE(false);
                }
            }
        }
    }
}

TEST_CASE("Grid: cubic_grid") {
    settings::grid::cubic = true;
    settings::grid::width = GENERATE(0.1, 0.5, 1);

    SECTION("largest x") {
        Limit3D axes(-10, 10, -1, 1, -1, 1);

        Grid grid(axes);
        auto gaxes = grid.get_axes();
        CHECK(gaxes.x.limits() == axes.x);
        CHECK(gaxes.x == gaxes.y);
        CHECK(gaxes.x == gaxes.z);
    }

    SECTION("largest y") {
        Limit3D axes(-1, 1, -10, 10, -1, 1);

        Grid grid(axes);
        auto gaxes = grid.get_axes();
        CHECK(gaxes.y.limits() == axes.y);
        CHECK(gaxes.y == gaxes.y);
        CHECK(gaxes.y == gaxes.z);
    }

    SECTION("largest x") {
        Limit3D axes(-1, 1, -1, 1, -10, 10);

        Grid grid(axes);
        auto gaxes = grid.get_axes();
        CHECK(gaxes.z.limits() == axes.z);
        CHECK(gaxes.z == gaxes.y);
        CHECK(gaxes.z == gaxes.z);
    }

    settings::grid::cubic = false;
}

TEST_CASE("Grid::operator=", "[files]") {
    settings::general::verbose = false;

    SECTION("copy") {
        Limit3D axes(-100, 100, -100, 100, -100, 100);
        Grid grid1(axes);
        {
            data::Molecule protein("test/files/2epe.pdb");
            grid1.add(protein.get_atoms());
            Grid grid2 = grid1;
            protein.set_grid(std::move(grid2));
            protein.generate_new_hydration();
            grid1.add(protein.get_waters());
        }

        Grid grid2 = grid1;
        REQUIRE(grid2 == grid1);

        Grid grid3(axes);
        grid3 = grid1;
        REQUIRE(grid3 == grid1);
    }

    SECTION("move") {
        Limit3D axes(-100, 100, -100, 100, -100, 100);
        Grid grid1(axes);
        {
            data::Molecule protein("test/files/2epe.pdb");
            grid1.add(protein.get_atoms());
            Grid grid2 = grid1;
            protein.set_grid(std::move(grid2));
            protein.generate_new_hydration();
            grid1.add(protein.get_waters());
        }

        Grid gridcopy = grid1;
        Grid grid2 = std::move(grid1);
        REQUIRE(grid2 == gridcopy);

        Grid grid3(axes);
        grid3 = std::move(grid2);
        REQUIRE(grid3 == gridcopy);
    }
}

TEST_CASE("Grid: hydration") {
    settings::general::verbose = false;

    Limit3D lims(-10, 10, -10, 10, -10, 10);
    vector<Atom> a = {Atom({0, 0, 0}, 0,  constants::atom_t::C, "LYS", 0)};
    data::Molecule protein(a);
    {
        Grid grid(lims);
        grid.add(a);
        protein.set_grid(std::move(grid));
        protein.generate_new_hydration();
    }

    // check that the waters have been expanded
    GridDebug* g = static_cast<GridDebug*>(protein.get_grid());
    for (const auto& w : g->get_member_waters()) {
        CHECK(w.is_expanded());
    }

    g->clear_waters();
    auto& gobj = g->grid;
    auto axes = g->get_axes();
    for (unsigned int i = 0; i < axes.x.bins; i++) {
        for (unsigned int j = 0; j < axes.y.bins; j++) {
            for (unsigned int k = 0; k < axes.z.bins; k++) {
                auto index = gobj.index(i, j, k);
                if (index != EMPTY) {
                    if (index == A_CENTER || index == A_AREA || index == VOLUME) {
                        continue;
                    }
                    std::cout << "Failed on index " << i << ", " << j << ", " << k << " with " << index << std::endl;
                    REQUIRE(false);
                }
            }
        }
    }
}

TEST_CASE("Grid::bin_ops") {
    Limit3D axes(-10, 10, -10, 10, -10, 10);
    Grid grid(axes);
    auto& gref = grid.grid;

    {
        gref.index(0, 0, 0) = A_CENTER;
        CHECK(gref.is_atom_center(0, 0, 0) == true);
        CHECK(gref.is_atom_area(0, 0, 0) == false);
        CHECK(gref.is_atom_area_or_volume(0, 0, 0) == false);
        CHECK(gref.is_water_center(0, 0, 0) == false);
        CHECK(gref.is_water_area(0, 0, 0) == false);
        CHECK(gref.is_empty(0, 0, 0) == false);
        CHECK(gref.is_volume(0, 0, 0) == false);
        CHECK(gref.is_empty_or_volume(0, 0, 0) == false);
    }
    {
        gref.index(0, 0, 0) = A_AREA;
        CHECK(gref.is_atom_center(0, 0, 0) == false);
        CHECK(gref.is_atom_area(0, 0, 0) == true);
        CHECK(gref.is_atom_area_or_volume(0, 0, 0) == true);
        CHECK(gref.is_water_center(0, 0, 0) == false);
        CHECK(gref.is_water_area(0, 0, 0) == false);
        CHECK(gref.is_empty(0, 0, 0) == false);
        CHECK(gref.is_volume(0, 0, 0) == false);
        CHECK(gref.is_empty_or_volume(0, 0, 0) == false);
    }
    {
        gref.index(0, 0, 0) = W_CENTER;
        CHECK(gref.is_atom_center(0, 0, 0) == false);
        CHECK(gref.is_atom_area(0, 0, 0) == false);
        CHECK(gref.is_atom_area_or_volume(0, 0, 0) == false);
        CHECK(gref.is_water_center(0, 0, 0) == true);
        CHECK(gref.is_water_area(0, 0, 0) == false);
        CHECK(gref.is_empty(0, 0, 0) == false);
        CHECK(gref.is_volume(0, 0, 0) == false);
        CHECK(gref.is_empty_or_volume(0, 0, 0) == false);
    }
    {
        gref.index(0, 0, 0) = W_AREA;
        CHECK(gref.is_atom_center(0, 0, 0) == false);
        CHECK(gref.is_atom_area(0, 0, 0) == false);
        CHECK(gref.is_atom_area_or_volume(0, 0, 0) == false);
        CHECK(gref.is_water_center(0, 0, 0) == false);
        CHECK(gref.is_water_area(0, 0, 0) == true);
        CHECK(gref.is_empty(0, 0, 0) == false);
        CHECK(gref.is_volume(0, 0, 0) == false);
        CHECK(gref.is_empty_or_volume(0, 0, 0) == false);
    }
    {
        gref.index(0, 0, 0) = VOLUME;
        CHECK(gref.is_atom_center(0, 0, 0) == false);
        CHECK(gref.is_atom_area(0, 0, 0) == false);
        CHECK(gref.is_atom_area_or_volume(0, 0, 0) == true);
        CHECK(gref.is_water_center(0, 0, 0) == false);
        CHECK(gref.is_water_area(0, 0, 0) == false);
        CHECK(gref.is_empty(0, 0, 0) == false);
        CHECK(gref.is_volume(0, 0, 0) == true);
        CHECK(gref.is_empty_or_volume(0, 0, 0) == true);
    }
    {
        gref.index(0, 0, 0) = EMPTY;
        CHECK(gref.is_atom_center(0, 0, 0) == false);
        CHECK(gref.is_atom_area(0, 0, 0) == false);
        CHECK(gref.is_atom_area_or_volume(0, 0, 0) == false);
        CHECK(gref.is_water_center(0, 0, 0) == false);
        CHECK(gref.is_water_area(0, 0, 0) == false);
        CHECK(gref.is_empty(0, 0, 0) == true);
        CHECK(gref.is_volume(0, 0, 0) == false);
        CHECK(gref.is_empty_or_volume(0, 0, 0) == true);
    }
}