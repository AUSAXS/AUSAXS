#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <grid/Grid.h>
#include <grid/detail/GridMember.h>
#include <hydrate/generation/RadialHydration.h>
#include <hydrate/generation/HydrationStrategy.h>
#include <hydrate/generation/HydrationFactory.h>
#include <hydrate/culling/CullingStrategy.h>
#include <settings/All.h>
#include <math/Vector3.h>
#include <settings/MoleculeSettings.h>
#include <constants/Constants.h>
#include <rigidbody/BodySplitter.h>

#include <vector>
#include <string> 

using namespace ausaxs;
using namespace ausaxs::grid;
using namespace ausaxs::grid::detail;
using namespace ausaxs::data;

// Debug class to expose the volume variable
class GridDebug : public grid::Grid {
    public: 
        ~GridDebug() override = default;
        GridDebug(Limit3D axes) : Grid(axes) {}

		double get_atomic_radius(form_factor::form_factor_t) const override {return ra;}
		double get_hydration_radius() const override {return rh;}
        void set_atomic_radius(double ra) {this->ra = ra;}
        void set_hydration_radius(double rh) {this->rh = rh;}
        auto get_volume_without_expanding() {return volume;}
        auto bounding_box_index() {return Grid::bounding_box_index();}
        auto to_bins(const Vector3<double> v) {return Grid::to_bins(v);}
        auto get_member_atoms() {return a_members;}
        auto get_member_waters() {return w_members;}
        static auto bounding_box(const std::vector<AtomFF>& atoms) {return Grid::bounding_box(atoms);}

    private:
        double ra = 0, rh = 0;
};

TEST_CASE("Grid::Grid") {
    settings::grid::cell_width = 1;

    SECTION("Limit3D&") {
        Limit3D axes(-10, 10, -10, 10, -10, 10);
        Grid grid(axes);
        CHECK(grid.a_members.empty());
        CHECK(grid.w_members.empty());
        CHECK(grid.a_members.empty());
        CHECK(grid.get_volume() == 0);
        CHECK(grid.get_width() == settings::grid::cell_width);
        CHECK(grid.get_axes() == Axis3D(axes, settings::grid::cell_width));
    }

    SECTION("vector<Atom>&") {
        std::vector<AtomFF> atoms = {
            AtomFF({0, 0, 0}, form_factor::form_factor_t::C), 
            AtomFF({1, 1, 1}, form_factor::form_factor_t::C), 
            AtomFF({2, 2, 2}, form_factor::form_factor_t::C)
        };
        Grid grid(atoms);
        CHECK(grid.a_members.size() == 3);
        CHECK(grid.w_members.empty());
        CHECK(grid.get_atoms() == atoms);
        CHECK(grid.get_volume() != 0);
        CHECK(grid.get_width() == settings::grid::cell_width);
    }

    SECTION("std::vector<Body>&") {
        std::vector<Body> bodies = {Body(std::vector{
            AtomFF({0, 0, 0}, form_factor::form_factor_t::C), 
            AtomFF({1, 1, 1}, form_factor::form_factor_t::C), 
            AtomFF({2, 2, 2}, form_factor::form_factor_t::C)
        })};
        Grid grid(bodies);
        CHECK(grid.a_members.size() == 3);
        CHECK(grid.w_members.empty());
        CHECK(grid.get_atoms() == bodies[0].get_atoms());
        CHECK(grid.get_volume() != 0);
        CHECK(grid.get_width() == settings::grid::cell_width);
    }

    SECTION("space_saving_constructor") {
        std::vector<AtomFF> atoms = {
            AtomFF({5,  0, -7}, form_factor::form_factor_t::C), 
            AtomFF({0, -5,  0}, form_factor::form_factor_t::C), 
            AtomFF({1,  1,  1}, form_factor::form_factor_t::C)
        };

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
            CHECK(axes.x.min == std::floor( 0 - 5*settings::grid::scaling*0.5 - settings::grid::cell_width));
            CHECK(axes.y.min == std::floor(-5 - 6*settings::grid::scaling*0.5 - settings::grid::cell_width));
            CHECK(axes.z.min == std::floor(-7 - 8*settings::grid::scaling*0.5 - settings::grid::cell_width));
            CHECK(axes.x.max ==  std::ceil(5 + 5*settings::grid::scaling*0.5 + settings::grid::cell_width));
            CHECK(axes.y.max ==  std::ceil(1 + 6*settings::grid::scaling*0.5 + settings::grid::cell_width));
            CHECK(axes.z.max ==  std::ceil(1 + 8*settings::grid::scaling*0.5 + settings::grid::cell_width));

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
            std::vector<Body> bodies; 
            for (const auto& a : atoms) {
                bodies.push_back(Body{std::vector{a}});
            }

            Grid grid(bodies);
            func(grid);
        }
    }
}

TEST_CASE("Grid::add") {
    Limit3D axes(-10, 10, -10, 10, -10, 10);
    settings::grid::cell_width = 1;
    GridDebug grid(axes);

    SECTION("Body") {
        SECTION("") {
            Body body{std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)}};
            grid.add(body, false);
            GridObj& g = grid.grid;

            CHECK(g.index(10, 10, 10) == A_CENTER);
            CHECK(g.index(10, 10, 11) == EMPTY);
            CHECK(g.index(10, 11, 11) == EMPTY);
            CHECK(g.index(11, 10, 10) == EMPTY);
            CHECK(g.index(9, 8, 14) == EMPTY);
            CHECK(grid.get_volume_without_expanding() == 1);
        }

        SECTION("") {
            Body body{std::vector{
                AtomFF({0, 0, 0}, form_factor::form_factor_t::C),
                AtomFF({1, 1, 1}, form_factor::form_factor_t::C)
            }};
            grid.add(body, false);
            GridObj& g = grid.grid;

            CHECK(g.index(10, 10, 10) == A_CENTER);
            CHECK(g.index(10, 10, 11) == EMPTY);
            CHECK(g.index(10, 11, 11) == EMPTY);
            CHECK(g.index(11, 10, 10) == EMPTY);
            CHECK(g.index(9, 8, 14) == EMPTY);

            CHECK(g.index(11, 11, 11) == A_CENTER);
            CHECK(g.index(11, 11, 12) == EMPTY);
            CHECK(g.index(11, 12, 12) == EMPTY);
            CHECK(g.index(12, 11, 11) == EMPTY);
            CHECK(g.index(10, 9, 15) == EMPTY);
            CHECK(grid.get_volume_without_expanding() == 2);
        }
    }

    // check that the grid is symmetry-aware
    SECTION("Body symmetries") {
        SECTION("") {
            Body body{std::vector{
                AtomFF({0, 0, 0}, form_factor::form_factor_t::C),
            }};
            body.symmetry().add(symmetry::Symmetry({{1, 1, 1}, {0, 0, 0}}));
            grid.add(body, false);
            GridObj& g = grid.grid;

            CHECK(g.index(10, 10, 10) == A_CENTER);
            CHECK(g.index(10, 10, 11) == EMPTY);
            CHECK(g.index(10, 11, 11) == EMPTY);
            CHECK(g.index(11, 10, 10) == EMPTY);
            CHECK(g.index(9, 8, 14) == EMPTY);

            CHECK(g.index(11, 11, 11) == A_CENTER);
            CHECK(g.index(11, 11, 12) == EMPTY);
            CHECK(g.index(11, 12, 12) == EMPTY);
            CHECK(g.index(12, 11, 11) == EMPTY);
            CHECK(g.index(10, 9, 15) == EMPTY);
            CHECK(grid.get_volume_without_expanding() == 2);
        }

        SECTION("") {
            Body body{std::vector{
                AtomFF({0, 0, 0}, form_factor::form_factor_t::C),
            }};
            body.symmetry().add(symmetry::Symmetry({{0, 0, 0}, {0, 0, 0}}, {{1, 1, 1}, {0, 0, 0}}, 5));
            grid.add(body, false);
            GridObj& g = grid.grid;

            CHECK(g.index(10, 10, 10) == A_CENTER);
            CHECK(g.index(11, 11, 11) == A_CENTER);
            CHECK(g.index(12, 12, 12) == A_CENTER);
            CHECK(g.index(13, 13, 13) == A_CENTER);
            CHECK(g.index(14, 14, 14) == A_CENTER);
            CHECK(g.index(15, 15, 15) == A_CENTER);
            CHECK(grid.get_volume_without_expanding() == 6);
        }
    }

    SECTION("Water") {
        SECTION("") {
            std::vector<Water> w = {Water({0, 0, 0})};
            for (auto& atom : w) {grid.add(atom, false);}
            GridObj& g = grid.grid;

            CHECK(g.index(10, 10, 10) == W_CENTER);
            CHECK(g.index(10, 10, 11) == EMPTY);
            CHECK(g.index(10, 11, 11) == EMPTY);
            CHECK(g.index(11, 10, 10) == EMPTY);
            CHECK(g.index(9, 8, 14) == EMPTY);
            CHECK(grid.get_volume() == 0);
        }

        SECTION("") {
            std::vector<Water> w = {
                Water({0, 0, 0}),
                Water({1, 1, 1})
            };
            for (auto& atom : w) {grid.add(atom, false);}
            GridObj& g = grid.grid;

            CHECK(g.index(10, 10, 10) == W_CENTER);
            CHECK(g.index(10, 10, 11) == EMPTY);
            CHECK(g.index(10, 11, 11) == EMPTY);
            CHECK(g.index(11, 10, 10) == EMPTY);
            CHECK(g.index(9, 8, 14) == EMPTY);

            CHECK(g.index(11, 11, 11) == W_CENTER);
            CHECK(g.index(11, 11, 12) == EMPTY);
            CHECK(g.index(11, 12, 12) == EMPTY);
            CHECK(g.index(12, 11, 11) == EMPTY);
            CHECK(g.index(10, 9, 15) == EMPTY);
            CHECK(grid.get_volume() == 0);
        }

        SECTION("") {
            std::vector<Water> w = {
                Water({0, 0, 0}),
                Water({1, 0, 0}),
                Water({0, 1, 0}),
                Water({0, 0, 1})
            };
            for (auto& atom : w) {grid.add(atom, false);}
            GridObj& g = grid.grid;

            CHECK(g.index(10, 10, 10) == W_CENTER);
            CHECK(g.index(11, 10, 10) == W_CENTER);
            CHECK(g.index(10, 11, 10) == W_CENTER);
            CHECK(g.index(10, 10, 11) == W_CENTER);
            CHECK(grid.get_volume() == 0);
        }
    }

    SECTION("std::vector<Water>") {
        SECTION("") {
            std::vector<Water> w = {Water({0, 0, 0})};
            grid.add(w, false);
            GridObj& g = grid.grid;

            CHECK(g.index(10, 10, 10) == W_CENTER);
            CHECK(g.index(10, 10, 11) == EMPTY);
            CHECK(g.index(10, 11, 11) == EMPTY);
            CHECK(g.index(11, 10, 10) == EMPTY);
            CHECK(g.index(9, 8, 14) == EMPTY);
            CHECK(grid.get_volume() == 0);
        }

        SECTION("") {
            std::vector<Water> w = {
                Water({0, 0, 0}),
                Water({1, 1, 1})
            };
            grid.add(w, false);
            GridObj& g = grid.grid;

            CHECK(g.index(10, 10, 10) == W_CENTER);
            CHECK(g.index(10, 10, 11) == EMPTY);
            CHECK(g.index(10, 11, 11) == EMPTY);
            CHECK(g.index(11, 10, 10) == EMPTY);
            CHECK(g.index(9, 8, 14) == EMPTY);

            CHECK(g.index(11, 11, 11) == W_CENTER);
            CHECK(g.index(11, 11, 12) == EMPTY);
            CHECK(g.index(11, 12, 12) == EMPTY);
            CHECK(g.index(12, 11, 11) == EMPTY);
            CHECK(g.index(10, 9, 15) == EMPTY);
            CHECK(grid.get_volume() == 0);
        }

        SECTION("") {
            std::vector<Water> w = {
                Water({0, 0, 0}),
                Water({1, 0, 0}),
                Water({0, 1, 0}),
                Water({0, 0, 1})
            };
            grid.add(w, false);
            GridObj& g = grid.grid;

            CHECK(g.index(10, 10, 10) == W_CENTER);
            CHECK(g.index(11, 10, 10) == W_CENTER);
            CHECK(g.index(10, 11, 10) == W_CENTER);
            CHECK(g.index(10, 10, 11) == W_CENTER);
            CHECK(grid.get_volume() == 0);
        }        
    }
}

TEST_CASE("Grid::expand_volume") {
    Limit3D lims(-10, 10, -10, 10, -10, 10);
    Grid grid(lims);

    SECTION("verify shape") {
        settings::grid::min_exv_radius = GENERATE(0.5, 1, 2);
        SECTION("no exv") {
            SECTION("") {
                settings::grid::min_exv_radius = 0;

                Body body{std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)}};
                grid.add(body, true);

                auto axes = grid.get_axes();
                double r = constants::radius::get_vdw_radius(form_factor::form_factor_t::C)/settings::grid::cell_width;
                for (unsigned int i = 0; i < axes.x.bins; ++i) {
                    for (unsigned int j = 0; j < axes.y.bins; ++j) {
                        for (unsigned int k = 0; k < axes.z.bins; ++k) {
                            if (i == 10 && j == 10 && k == 10) {continue;}
                            double dist = grid.to_xyz(i, j, k).norm(); // dist from center
                            if (dist <= r) {
                                REQUIRE(grid.grid.index(i, j, k) == A_AREA);
                            } else {
                                REQUIRE(grid.grid.index(i, j, k) == EMPTY);
                            }
                        }
                    }
                }
                CHECK(grid.grid.index(10, 10, 10) == A_CENTER);
            }

            SECTION("") {
                settings::grid::min_exv_radius = 0;

                Body body{std::vector{
                    AtomFF({ 2, 0, 0}, form_factor::form_factor_t::C),
                    AtomFF({-2, 0, 0}, form_factor::form_factor_t::N)
                }};
                grid.add(body, true);

                auto axes = grid.get_axes();
                double rC = constants::radius::get_vdw_radius(form_factor::form_factor_t::C)/settings::grid::cell_width;
                double rN = constants::radius::get_vdw_radius(form_factor::form_factor_t::N)/settings::grid::cell_width;
                for (unsigned int i = 0; i < axes.x.bins; ++i) {
                    for (unsigned int j = 0; j < axes.y.bins; ++j) {
                        for (unsigned int k = 0; k < axes.z.bins; ++k) {
                            double dist1 = grid.to_xyz(i, j, k).distance(body.get_atom(0).coordinates());
                            double dist2 = grid.to_xyz(i, j, k).distance(body.get_atom(1).coordinates());
                            if (dist1 == 0 || dist2 == 0) {continue;}
                            if (dist1 <= rC || dist2 <= rN) {
                                REQUIRE(grid.grid.index(i, j, k) == A_AREA);
                            } else {
                                REQUIRE(grid.grid.index(i, j, k) == EMPTY);
                            }
                        }
                    }
                }
                CHECK(grid.grid.index(grid.to_bins(body.get_atom(0).coordinates())) == A_CENTER);
                CHECK(grid.grid.index(grid.to_bins(body.get_atom(1).coordinates())) == A_CENTER);
            }
        }

        SECTION("exv") {
            settings::grid::min_exv_radius = 3;

            Body body{std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)}};
            grid.add(body, true);

            auto axes = grid.get_axes();
            double r = constants::radius::get_vdw_radius(form_factor::form_factor_t::C)/settings::grid::cell_width;
            for (unsigned int i = 0; i < axes.x.bins; ++i) {
                for (unsigned int j = 0; j < axes.y.bins; ++j) {
                    for (unsigned int k = 0; k < axes.z.bins; ++k) {
                        if (i == 10 && j == 10 && k == 10) {continue;}
                        double dist = grid.to_xyz(i, j, k).norm(); // dist from center
                        if (dist <= r) {
                            REQUIRE(grid.grid.index(i, j, k) == A_AREA);
                        } else if (dist <= 3) {
                            REQUIRE(grid.grid.index(i, j, k) == VOLUME);
                        } else {
                            REQUIRE(grid.grid.index(i, j, k) == EMPTY);
                        }
                    }
                }
            }
            CHECK(grid.grid.index(10, 10, 10) == A_CENTER);
        }
    }

    // these atom counts have been checked by visual inspection
    SECTION("count tests") {
        settings::grid::min_bins = 30;

        SECTION("five atoms") {
            settings::grid::min_exv_radius = 3;
            Body body{std::vector{
                AtomFF({ 0,  0,  0}, form_factor::form_factor_t::C),
                AtomFF({-2,  0,  0}, form_factor::form_factor_t::C),
                AtomFF({ 2,  0,  0}, form_factor::form_factor_t::C),
                AtomFF({ 0, -2,  0}, form_factor::form_factor_t::C),
                AtomFF({ 0,  2,  0}, form_factor::form_factor_t::C)
            }};

            GridDebug grid({-10, 10, -10, 10, -10, 10});
            grid.add(body, true);
            REQUIRE(grid.get_volume() == 323);
        }
        settings::grid::min_bins = 0;
    }
}

TEST_CASE("Grid::remove") {
    Limit3D axes(-10, 10, -10, 10, -10, 10);
    settings::grid::cell_width = 1;
    Grid grid(axes);

    SECTION("") {
        Body b{std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)}};
        grid.add(b, true);
        grid.remove(b);
        GridObj& g = grid.grid;

        auto axes = grid.get_axes();
        for (unsigned int i = 0; i < axes.x.bins; i++) {
            for (unsigned int j = 0; j < axes.y.bins; j++) {
                for (unsigned int k = 0; k < axes.z.bins; k++) {
                    REQUIRE(g.index(i, j, k) == EMPTY);
                }
            }
        }
        REQUIRE(grid.get_volume() == 0);
    }

    SECTION("") {
        Body b{std::vector{
            AtomFF({0, 0, 0}, form_factor::form_factor_t::C),
            AtomFF({1, 1, 1}, form_factor::form_factor_t::C)
        }};
        grid.add(b, true);
        grid.remove(b);
        GridObj& g = grid.grid;

        auto axes = grid.get_axes();
        for (unsigned int i = 0; i < axes.x.bins; i++) {
            for (unsigned int j = 0; j < axes.y.bins; j++) {
                for (unsigned int k = 0; k < axes.z.bins; k++) {
                    REQUIRE(g.index(i, j, k) == EMPTY);
                }
            }
        }
        REQUIRE(grid.get_volume() == 0);
    }

    SECTION("symmetry-aware") {
        Body b{std::vector{
            AtomFF({0, 0, 0}, form_factor::form_factor_t::C),
        }};
        b.symmetry().add(symmetry::Symmetry({{1, 1, 1}, {0, 0, 0}}));
        grid.add(b, true);
        grid.remove(b);
        GridObj& g = grid.grid;

        auto axes = grid.get_axes();
        for (unsigned int i = 0; i < axes.x.bins; i++) {
            for (unsigned int j = 0; j < axes.y.bins; j++) {
                for (unsigned int k = 0; k < axes.z.bins; k++) {
                    REQUIRE(g.index(i, j, k) == EMPTY);
                }
            }
        }
        REQUIRE(grid.get_volume() == 0);
    }
}

TEST_CASE("Grid::remove_waters") {
    Limit3D axes(-10, 10, -10, 10, -10, 10);
    settings::grid::cell_width = 1;
    Grid grid(axes);

    SECTION("") {
        Water w({0, 0, 0});
        grid.add(w, true);
        grid.remove_waters({true});
        GridObj& g = grid.grid;

        REQUIRE(grid.get_waters().empty());
        auto axes = grid.get_axes();
        for (unsigned int i = 0; i < axes.x.bins; i++) {
            for (unsigned int j = 0; j < axes.y.bins; j++) {
                for (unsigned int k = 0; k < axes.z.bins; k++) {
                    REQUIRE(g.index(i, j, k) == EMPTY);
                }
            }
        }
    }

    SECTION("") {
        auto w = {
            Water({0, 0, 0}),
            Water({1, 0, 0}),
            Water({0, 1, 0}),
            Water({0, 0, 1})
        };
        grid.add(w, true);
        grid.remove_waters({true, true, true, true});
        GridObj& g = grid.grid;

        REQUIRE(grid.get_waters().empty());
        auto axes = grid.get_axes();
        for (unsigned int i = 0; i < axes.x.bins; i++) {
            for (unsigned int j = 0; j < axes.y.bins; j++) {
                for (unsigned int k = 0; k < axes.z.bins; k++) {
                    REQUIRE(g.index(i, j, k) == EMPTY);
                }
            }
        }
    }

    SECTION("") {
        std::vector w = {
            Water({0, 0, 0}),
            Water({1, 0, 0}),
            Water({0, 1, 0}),
            Water({0, 0, 1})
        };
        grid.add(w, true);
        grid.remove_waters({true, true, true, false});

        auto wg = grid.get_waters();
        REQUIRE(wg.size() == 1);
        REQUIRE(wg[0] == w[3]);
    }
}

TEST_CASE("Grid: GridMember") {
    grid::GridMember<AtomFF> a(AtomFF({0.1, 0.2, 0.3}, form_factor::form_factor_t::C), {1, 2, 3});

    auto b = a;
    CHECK(b == a);

    auto c(std::move(b));
    CHECK(c == a);
}

TEST_CASE("Grid::bounding_box") {
    Limit3D axes(-10, 10, -10, 10, -10, 10);
    settings::grid::cell_width = 1;
    GridDebug grid(axes);

    SECTION("simple") {
        Body body{std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)}};
        grid.add(body, false);

        auto[min, max] = grid.bounding_box_index();
        CHECK(min[0] == 10);
        CHECK(max[0] == 11);
        CHECK(min[1] == 10);
        CHECK(max[1] == 11);
        CHECK(min[2] == 10);
        CHECK(max[2] == 11);
    }

    SECTION("complex") {
        Body body{std::vector{
            AtomFF({5, 0, -7}, form_factor::form_factor_t::C), 
            AtomFF({0, -5, 0}, form_factor::form_factor_t::C)
        }};
        grid.add(body);
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
    SECTION("simple") {
        Limit3D lims(-10, 10, -10, 10, -10, 10);
        settings::grid::cell_width = 1e-1;
        Grid grid(lims);

        Body body{std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)}};
        grid.add(body, true);
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
        CHECK(grid.get_volume() == count*std::pow(settings::grid::cell_width, 3));
    }

    SECTION("multiple atoms") {
        Limit3D axes(-10, 10, -10, 10, -10, 10);
        settings::grid::cell_width = 1;
        GridDebug grid(axes);
        grid.set_atomic_radius(1);

        // cout << grid.get_volume() << endl;
        Body body{std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)}};
        std::vector<Water> w = {
            Water({2, 2, 2}), 
            Water({2, 2, 3})
        };
        grid.add(body, false);
        grid.add(w, false);
        REQUIRE(grid.get_volume_without_expanding() == 1); // atoms are added as point-particles, and only occupy one unit of space.

        SECTION("") {
            settings::grid::min_exv_radius = 1;
            REQUIRE(grid.get_volume() == 7); // the radius is 1, so expanding the volume in a sphere results in one unit of volume added along each coordinate axis

            grid.add(Body{std::vector{{AtomFF({0, 0, -1}, form_factor::form_factor_t::C)}}});
            REQUIRE(grid.get_volume() == 12); // second atom is placed adjacent to the first one, so the volumes overlap. 
        }

        SECTION("") {
            settings::grid::min_exv_radius = 2;
            grid.expand_volume();
            REQUIRE(grid.get_volume() == 33);
        }
    }
}

TEST_CASE("Grid::width") {
    // check that the grid width is actually used
    settings::grid::cell_width = GENERATE(0.25, 0.5, 1, 2);

    Limit3D lims(-10, 10, -10, 10, -10, 10);
    Grid grid(lims);

    auto axes = grid.get_axes();
    for (unsigned int i = 0; i < axes.x.bins-1; ++i) {
        CHECK(grid.to_xyz(i, 0, 0).distance(grid.to_xyz(i+1, 0, 0)) == settings::grid::cell_width);
        CHECK(grid.to_xyz(0, i, 0).distance(grid.to_xyz(0, i+1, 0)) == settings::grid::cell_width);
        CHECK(grid.to_xyz(0, 0, i).distance(grid.to_xyz(0, 0, i+1)) == settings::grid::cell_width);
    }

    auto center = grid.get_center();
    CHECK(grid.to_xyz(center) == Vector3<double>{0, 0, 0});
    CHECK(grid.to_xyz(center+Vector3<double>{1, 0, 0}) == Vector3<double>{settings::grid::cell_width, 0, 0});
    CHECK(grid.to_xyz(center+Vector3<double>{0, 1, 0}) == Vector3<double>{0, settings::grid::cell_width, 0});
    CHECK(grid.to_xyz(center+Vector3<double>{0, 0, 1}) == Vector3<double>{0, 0, settings::grid::cell_width});
}

TEST_CASE("Grid::hydrate") {
    settings::general::verbose = false;
    settings::molecule::center = false;
    settings::hydrate::shell_correction = 0;
    settings::hydrate::hydration_strategy = GENERATE(
        settings::hydrate::HydrationStrategy::AxesStrategy, 
        settings::hydrate::HydrationStrategy::RadialStrategy,
        settings::hydrate::HydrationStrategy::JanStrategy
    );
    hydrate::RadialHydration::set_noise_generator([] () {return Vector3<double>{0, 0, 0};});

    // check that all the expected hydration sites are found
    SECTION("correct placement " + std::to_string(static_cast<int>(settings::hydrate::hydration_strategy))) {
        Limit3D axes(-10, 10, -10, 10, -10, 10);
        settings::grid::cell_width = 1;
        auto grid = std::make_unique<GridDebug>(axes);
        grid->set_atomic_radius(3);
        grid->set_hydration_radius(3);

        // add a single atom to the grid, and hydrate it
        settings::grid::water_scaling = 0;
        Body b{std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)}};
        grid->add(b);
        REQUIRE(grid->get_atoms().size() == 1); // check atoms

        data::Molecule protein({b});
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
    SECTION("reversible " + std::to_string(static_cast<int>(settings::hydrate::hydration_strategy))) {
        SECTION("LAR1-2") {
            // get the grid before hydrating it
            Molecule protein("tests/files/LAR1-2.pdb");
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
        Molecule protein("tests/files/LAR1-2.pdb");
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
    settings::hydrate::shell_correction = 0;
    hydrate::RadialHydration::set_noise_generator([] () {return Vector3<double>{0, 0, 0};});

    auto test_width_basics = [] (settings::hydrate::HydrationStrategy strategy) {
        settings::hydrate::hydration_strategy = strategy;
        settings::grid::cell_width = GENERATE(
            0.5,
            1
        );
        settings::grid::min_exv_radius = 0;
        Body b{std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)}};
        data::Molecule protein({b});
        {
            Limit3D axes(-10, 10, -10, 10, -10, 10);
            auto grid = std::make_unique<GridDebug>(axes);
            grid->set_atomic_radius(3);
            grid->set_hydration_radius(3);

            grid->add(b, false);
            GridObj& g = grid->grid;

            // check that it was placed correctly
            double icw = 1./settings::grid::cell_width;
            REQUIRE(g.index(10*icw, 10*icw, 10*icw) == A_CENTER);
            REQUIRE(grid->a_members.back().get_bin_loc() == Vector3(10*icw, 10*icw, 10*icw));
            REQUIRE(grid->a_members.back().get_atom().coordinates() == Vector3(0, 0, 0));

            // check water generation
            settings::grid::water_scaling = 0;

            protein.set_grid(std::move(grid));
            protein.generate_new_hydration();
        }
        std::vector<Water> water = protein.get_waters();

        // since this needs to work with different placement strategies, we have to perform a more general check on the positions
        std::vector<Vector3<int>> v = {{0, 0, 6}, {0, 0, -6}, {6, 0, 0}, {-6, 0, 0}, {0, 6, 0}, {0, -6, 0}};
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
        settings::grid::cell_width = 0.1;
        Limit3D axes(-10, 10, -10, 10, -10, 10);
        GridDebug grid(axes);
        grid.set_atomic_radius(3);
        grid.set_hydration_radius(3);

        Body b{std::vector{
            AtomFF({5, 0, -7}, form_factor::form_factor_t::C), 
            AtomFF({0, -5, 0}, form_factor::form_factor_t::C)
        }};
        grid.add(b, true);

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
    settings::grid::cell_width = 1;
    GridDebug grid(axes);
    grid.set_atomic_radius(3);
    grid.set_hydration_radius(3);

    // atoms
    AtomFF a1({3, 0, 0}, form_factor::form_factor_t::C);
    AtomFF a2({0, 3, 0}, form_factor::form_factor_t::C);
    AtomFF a3({0, 0, 3}, form_factor::form_factor_t::C);

    // waters
    Water w1(Vector3<double>{0, 0, -3});
    Water w2(Vector3<double>{0, -3, 0});
    Water w3(Vector3<double>{-3, 0, 0});
    Body body{std::vector{a1, a2, a3}, std::vector{w1, w2, w3}};

    SECTION("add") {
        // add atoms
        grid.add(body);
        auto ga = grid.get_atoms();
        REQUIRE(grid.a_members.size() == 3); // check actual data
        REQUIRE(ga.size() >= 3);
        auto a = body.get_atoms();
        for (int i = 0; i < 3; i++) {
            REQUIRE(ga[i] == a[i]); // check equality with what we added
        }

        // add waters
        std::vector<Water> wa = grid.get_waters();
        REQUIRE(grid.a_members.size() == 3); // check actual data
        REQUIRE(grid.w_members.size() == 3); // check actual data
        REQUIRE(wa.size() >= 3);
        auto w = body.get_waters()->get();
        for (int i = 0; i < 3; i++) {
            REQUIRE(wa[i] == w[i]); // check equality with what we added
        }
    }

    SECTION("remove") {
        grid.add(body);
        grid.remove(body);
        auto ga = grid.get_atoms();
        auto wa = grid.get_waters();

        // check sizes
        REQUIRE(ga.size() == 0);
        REQUIRE(wa.size() == 0);

        // check old centers
        GridObj &g = grid.grid;
        for (auto a : body.get_atoms()) {
            auto loc = grid.to_bins(a.coordinates());
            REQUIRE(g.index(loc) == EMPTY);
        }
    }

    SECTION("remove index") {
        SECTION("water") {
            grid.add(body);

            REQUIRE(grid.w_members.size() == 3);
            grid.remove_waters({true, true, false});
            REQUIRE(grid.w_members.size() == 1);

            // check the remaining one
            auto ga = grid.get_atoms();
            auto wa = grid.get_waters();
            REQUIRE(wa.size() == 1);
            REQUIRE(ga.size() == 3);
            REQUIRE(wa[0] == w3);
        }
    }

    SECTION("clear_waters") {
        grid.add(body);
        grid.clear_waters();
        REQUIRE(grid.a_members.size() == 3);
        REQUIRE(grid.w_members.size() == 0);
    }
}

TEST_CASE("Grid: correct_volume") {
    Limit3D axes(-10, 10, -10, 10, -10, 10);
    settings::grid::cell_width = 1;
    GridDebug grid(axes);
    grid.set_atomic_radius(10);     // heavy overlap
    grid.set_hydration_radius(10);

    // atoms
    AtomFF a1({3, 0, 0}, form_factor::form_factor_t::C);
    AtomFF a2({0, 3, 0}, form_factor::form_factor_t::C);
    AtomFF a3({0, 0, 3}, form_factor::form_factor_t::C);
    std::vector<Body> body{
        Body{std::vector{{a1}}}, 
        Body{std::vector{{a2}}}, 
        Body{std::vector{{a3}}}
    };

    // waters
    Water w1({0, 0, -3});
    Water w2({0, -3, 0});
    Water w3({-3, 0, 0});
    std::vector<Body> w = {
        {{}, std::vector{w1}}, 
        {{}, std::vector{w2}}, 
        {{}, std::vector{w3}}
    };

    // single non-overlapping
    REQUIRE(grid.get_volume_without_expanding() == 0);
    grid.add(body[0]);
    REQUIRE(grid.get_volume_without_expanding() != 0);
    grid.remove(body[0]);
    REQUIRE(grid.get_volume_without_expanding() == 0);

    grid.add(w[0]);
    REQUIRE(grid.get_volume_without_expanding() == 0);
    grid.remove(w[0]);
    REQUIRE(grid.get_volume_without_expanding() == 0);

    // multiple overlapping
    for(auto& b : body) {grid.add(b);}
    for (auto& b : w) {grid.add(b);}
    REQUIRE(grid.get_volume_without_expanding() != 0);
    for(auto& b : body) {grid.remove(b);}
    for (auto& b : w) {grid.remove(b);}
    REQUIRE(grid.get_volume_without_expanding() == 0);
}

TEST_CASE("Grid::find_free_locs") {
    settings::general::verbose = false;
    settings::hydrate::shell_correction = 0;
    hydrate::RadialHydration::set_noise_generator([] () {return Vector3<double>{0, 0, 0};});

    auto test_func = [] (settings::hydrate::HydrationStrategy ch) {
        settings::hydrate::hydration_strategy = ch;
        settings::grid::cell_width = 1;

        Limit3D axes(-10, 10, -10, 10, -10, 10);
        auto grid = std::make_unique<GridDebug>(axes);
        grid->set_atomic_radius(3);
        grid->set_hydration_radius(3);

        settings::grid::water_scaling = 0;
        Body body{std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)}};
        grid->add(body);
        grid->expand_volume();

        data::Molecule protein({body});
        protein.set_grid(std::move(grid));
        protein.generate_new_hydration();

        std::vector<grid::GridMember<Water>> locs = static_cast<GridDebug*>(protein.get_grid())->get_member_waters();
        REQUIRE(locs.size() == 6);

        // since this needs to work with different placement strategies, we have to perform a more general check on the positions
        std::vector<Vector3<int>> v = {{0, 0, 6}, {0, 0, -6}, {6, 0, 0}, {-6, 0, 0}, {0, 6, 0}, {0, -6, 0}};
        for (const auto& l : locs) {
            bool found = false;
            for (const auto& p : v) {
                if (l.get_atom().coordinates() == p) {found = true;}
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
    settings::grid::min_bins = 0;

    Limit3D axes(-10, 10, -10, 10, -10, 10);
    settings::grid::cell_width = 1;
    auto grid = std::make_unique<GridDebug>(axes);
    grid->set_atomic_radius(3);
    grid->set_hydration_radius(3);

    Body body{std::vector{
        AtomFF({3, 0, 0}, form_factor::form_factor_t::C), 
        AtomFF({0, 3, 0}, form_factor::form_factor_t::C)
    }};
    grid->add(body, false);
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
    settings::grid::cell_width = GENERATE(0.1, 0.5, 1);

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
            data::Molecule protein("tests/files/2epe.pdb");
            protein.clear_hydration();
            grid1.add(protein.get_body(0));
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
            data::Molecule protein("tests/files/2epe.pdb");
            protein.clear_hydration();
            grid1.add(protein.get_body(0));
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
    Body body{std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)}};
    data::Molecule protein({body});
    {
        Grid grid(lims);
        grid.add(body);
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

TEST_CASE("Grid::add:remove") {
    SECTION("single") {
        Body b(std::vector<AtomFF>{
            AtomFF({-1, -1, -1}, form_factor::form_factor_t::C), 
            AtomFF({-1,  1, -1}, form_factor::form_factor_t::C)
        });
        grid::Grid g(Limit3D(-2, 2, -2, 2, -2, 2));

        g.add(b);
        REQUIRE(g.a_members.size() == 2);
        REQUIRE(g.get_volume() != 0);

        g.remove(b);
        REQUIRE(g.a_members.size() == 0);
        REQUIRE(g.get_volume() == 0);
    }

    SECTION("multiple") {
        constants::radius::set_dummy_radius(1);
        std::vector<AtomFF> a1 = {AtomFF({-1, -1, -1}, form_factor::form_factor_t::C), AtomFF({-1, 1, -1}, form_factor::form_factor_t::C)};
        std::vector<AtomFF> a2 = {AtomFF({ 1, -1, -1}, form_factor::form_factor_t::C), AtomFF({ 1, 1, -1}, form_factor::form_factor_t::C)};
        std::vector<AtomFF> a3 = {AtomFF({-1, -1,  1}, form_factor::form_factor_t::C), AtomFF({-1, 1,  1}, form_factor::form_factor_t::C)};
        std::vector<AtomFF> a4 = {AtomFF({ 1, -1,  1}, form_factor::form_factor_t::C), AtomFF({ 1, 1,  1}, form_factor::form_factor_t::C)};
        Body b1(a1), b2(a2), b3(a3), b4(a4);
        std::vector<Body> bodies = {b1, b2, b3, b4};
        grid::Grid grid(Limit3D(-5, 5, -5, 5, -5, 5));

        grid.add(b1);
        grid.add(b2);
        grid.add(b3);
        grid.add(b4);
        REQUIRE(grid.a_members.size() == 8);

        unsigned int vol = grid.get_volume();
        grid.remove(b2);
        grid.add(b2);
        REQUIRE(grid.get_volume() == vol);
        grid.remove(b2);
        grid.force_expand_volume();
        REQUIRE(grid.a_members.size() == 6);

        vol = grid.get_volume();
        grid.remove(b1);
        grid.add(b1);
        REQUIRE(grid.get_volume() == vol);
        grid.remove(b1);
        grid.force_expand_volume();
        REQUIRE(grid.a_members.size() == 4);

        vol = grid.get_volume();
        grid.remove(b3);
        grid.add(b3);
        REQUIRE(grid.get_volume() == vol);
        grid.remove(b3);
        grid.force_expand_volume();
        REQUIRE(grid.a_members.size() == 2);

        auto remaining = grid.a_members;
        for (const auto& e : remaining) {
            REQUIRE((e == a4[0] || e == a4[1]));
        }

        // check volume
        REQUIRE(grid.get_volume() != 0);
        grid.remove(b4);
        REQUIRE(grid.get_volume() == 0);
    }

    SECTION("real data") {
        settings::general::verbose = false;
        Molecule protein = rigidbody::BodySplitter::split("tests/files/2epe.pdb", {9, 99});
        unsigned int N = protein.get_atoms().size();
        auto grid = protein.get_grid();
        CHECK(grid->a_members.size() == N);
        CHECK(grid->get_volume() != 0);

        // body 1
        int vol = grid->get_volume();
        grid->remove(protein.get_body(0));
        grid->add(   protein.get_body(0));
        CHECK(grid->get_volume() == vol);

        grid->remove(protein.get_body(0));
        grid->force_expand_volume();
        CHECK(grid->a_members.size() == N - protein.get_body(0).size_atom());
        CHECK(grid->get_volume() != 0);

        // body 2
        vol = grid->get_volume();
        grid->remove(protein.get_body(1));
        grid->add(   protein.get_body(1));
        CHECK(grid->get_volume() == vol);

        grid->remove(protein.get_body(1));        
        grid->force_expand_volume();
        CHECK(grid->a_members.size() == N - protein.get_body(0).size_atom() - protein.get_body(1).size_atom());
        CHECK(grid->get_volume() != 0);

        // body 3
        vol = grid->get_volume();
        grid->remove(protein.get_body(2));
        grid->add(   protein.get_body(2));
        CHECK(grid->get_volume() == vol);

        grid->remove(protein.get_body(2));
        CHECK(grid->a_members.size() == 0);
        CHECK(grid->get_volume() == 0);
    }
}