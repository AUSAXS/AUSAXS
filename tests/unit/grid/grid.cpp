#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <grid/Grid.h>
#include <grid/detail/GridMember.h>
#include <data/Body.h>
#include <data/atoms/AtomFF.h>
#include <data/atoms/Water.h>
#include <settings/GridSettings.h>
#include <utility/Axis3D.h>

using namespace ausaxs;
using namespace ausaxs::grid;
using namespace ausaxs::grid::detail;
using namespace ausaxs::data;

TEST_CASE("Grid::constructor Limit3D") {
    settings::grid::cell_width = 1.0;
    Limit3D axes(-10, 10, -10, 10, -10, 10);
    Grid grid(axes);

    SECTION("empty grid") {
        REQUIRE(grid.a_members.empty());
        REQUIRE(grid.w_members.empty());
        REQUIRE(grid.get_volume_bins() == 0);
    }

    SECTION("axes") {
        REQUIRE(grid.get_axes() == Axis3D(axes, settings::grid::cell_width));
        REQUIRE(grid.get_width() == settings::grid::cell_width);
    }

    SECTION("bins") {
        auto bins = grid.get_bins();
        REQUIRE(bins.x() == 20);
        REQUIRE(bins.y() == 20);
        REQUIRE(bins.z() == 20);
    }
}

TEST_CASE("Grid::constructor vector<AtomFF>") {
    settings::grid::cell_width = 1.0;
    std::vector<AtomFF> atoms = {
        AtomFF({0, 0, 0}, form_factor::form_factor_t::C),
        AtomFF({1, 1, 1}, form_factor::form_factor_t::C),
        AtomFF({2, 2, 2}, form_factor::form_factor_t::C)
    };
    Grid grid(atoms);

    SECTION("atoms added") {
        REQUIRE(grid.a_members.size() == 3);
        REQUIRE(grid.w_members.empty());
    }

    SECTION("get_atoms") {
        auto retrieved_atoms = grid.get_atoms();
        REQUIRE(retrieved_atoms.size() == 3);
        for (size_t i = 0; i < atoms.size(); ++i) {
            REQUIRE(retrieved_atoms[i] == atoms[i]);
        }
    }
}

TEST_CASE("Grid::constructor vector<Body>") {
    settings::grid::cell_width = 1.0;
    std::vector<AtomFF> atoms = {
        AtomFF({0, 0, 0}, form_factor::form_factor_t::C),
        AtomFF({1, 1, 1}, form_factor::form_factor_t::C)
    };
    std::vector<Body> bodies = {Body(atoms)};
    Grid grid(bodies);

    SECTION("atoms added") {
        REQUIRE(grid.a_members.size() == 2);
        REQUIRE(grid.w_members.empty());
    }

    SECTION("get_atoms") {
        auto retrieved_atoms = grid.get_atoms();
        REQUIRE(retrieved_atoms == atoms);
    }
}

TEST_CASE("Grid::copy constructor") {
    settings::grid::cell_width = 1.0;
    Limit3D axes(-5, 5, -5, 5, -5, 5);
    Grid grid1(axes);
    
    AtomFF atom({0, 0, 0}, form_factor::form_factor_t::C);
    Body body{std::vector{atom}};
    grid1.add(body, false);

    Grid grid2(grid1);

    REQUIRE(grid2.a_members.size() == grid1.a_members.size());
    REQUIRE(grid2.get_axes() == grid1.get_axes());
    REQUIRE(grid2.get_width() == grid1.get_width());
}

TEST_CASE("Grid::move constructor") {
    settings::grid::cell_width = 1.0;
    Limit3D axes(-5, 5, -5, 5, -5, 5);
    Grid grid1(axes);
    
    AtomFF atom({0, 0, 0}, form_factor::form_factor_t::C);
    Body body{std::vector{atom}};
    grid1.add(body, false);
    
    size_t original_size = grid1.a_members.size();
    Grid grid2(std::move(grid1));

    REQUIRE(grid2.a_members.size() == original_size);
}

TEST_CASE("Grid::operator=") {
    settings::grid::cell_width = 1.0;
    Limit3D axes1(-5, 5, -5, 5, -5, 5);
    Limit3D axes2(-10, 10, -10, 10, -10, 10);
    Grid grid1(axes1);
    Grid grid2(axes2);
    
    AtomFF atom({0, 0, 0}, form_factor::form_factor_t::C);
    Body body{std::vector{atom}};
    grid1.add(body, false);

    SECTION("copy assignment") {
        grid2 = grid1;
        REQUIRE(grid2.a_members.size() == grid1.a_members.size());
        REQUIRE(grid2.get_axes() == grid1.get_axes());
    }

    SECTION("move assignment") {
        size_t original_size = grid1.a_members.size();
        grid2 = std::move(grid1);
        REQUIRE(grid2.a_members.size() == original_size);
    }
}

TEST_CASE("Grid::operator==") {
    settings::grid::cell_width = 1.0;
    Limit3D axes(-5, 5, -5, 5, -5, 5);
    Grid grid1(axes);
    Grid grid2(axes);

    SECTION("equal grids") {
        REQUIRE(grid1 == grid2);
    }

    SECTION("different axes") {
        Limit3D different_axes(-10, 10, -10, 10, -10, 10);
        Grid grid3(different_axes);
        REQUIRE_FALSE(grid1 == grid3);
    }
}

TEST_CASE("Grid::get_bins") {
    settings::grid::cell_width = 2.0;
    Limit3D axes(-10, 10, -20, 20, -30, 30);
    Grid grid(axes);

    auto bins = grid.get_bins();
    REQUIRE(bins.x() == 10);
    REQUIRE(bins.y() == 20);
    REQUIRE(bins.z() == 30);
}

TEST_CASE("Grid::get_width") {
    settings::grid::cell_width = 1.5;
    Limit3D axes(-10, 10, -10, 10, -10, 10);
    Grid grid(axes);

    REQUIRE_THAT(grid.get_width(), Catch::Matchers::WithinAbs(1.5, 1e-6));
}

TEST_CASE("Grid::to_bins") {
    settings::grid::cell_width = 1.0;
    Limit3D axes(-10, 10, -10, 10, -10, 10);
    Grid grid(axes);

    SECTION("origin") {
        auto bins = grid.to_bins({0, 0, 0});
        REQUIRE(bins.x() == 10);
        REQUIRE(bins.y() == 10);
        REQUIRE(bins.z() == 10);
    }

    SECTION("positive coordinates") {
        auto bins = grid.to_bins({5, 5, 5});
        REQUIRE(bins.x() == 15);
        REQUIRE(bins.y() == 15);
        REQUIRE(bins.z() == 15);
    }

    SECTION("negative coordinates") {
        auto bins = grid.to_bins({-5, -5, -5});
        REQUIRE(bins.x() == 5);
        REQUIRE(bins.y() == 5);
        REQUIRE(bins.z() == 5);
    }
}

TEST_CASE("Grid::to_bins_bounded") {
    settings::grid::cell_width = 1.0;
    Limit3D axes(-10, 10, -10, 10, -10, 10);
    Grid grid(axes);

    SECTION("within bounds") {
        auto bins = grid.to_bins_bounded({0, 0, 0});
        REQUIRE(bins.x() == 10);
        REQUIRE(bins.y() == 10);
        REQUIRE(bins.z() == 10);
    }

    SECTION("outside bounds - clamped") {
        auto bins = grid.to_bins_bounded({100, -100, 50});
        REQUIRE(bins.x() >= 0);
        REQUIRE(bins.x() < grid.get_bins().x());
        REQUIRE(bins.y() >= 0);
        REQUIRE(bins.y() < grid.get_bins().y());
        REQUIRE(bins.z() >= 0);
        REQUIRE(bins.z() < grid.get_bins().z());
    }
}

TEST_CASE("Grid::to_xyz") {
    settings::grid::cell_width = 1.0;
    Limit3D axes(-10, 10, -10, 10, -10, 10);
    Grid grid(axes);

    SECTION("from Vector3") {
        auto xyz = grid.to_xyz(Vector3<int>(10, 10, 10));
        REQUIRE_THAT(xyz.x(), Catch::Matchers::WithinAbs(0.0, 1e-6));
        REQUIRE_THAT(xyz.y(), Catch::Matchers::WithinAbs(0.0, 1e-6));
        REQUIRE_THAT(xyz.z(), Catch::Matchers::WithinAbs(0.0, 1e-6));
    }

    SECTION("from three integers") {
        auto xyz = grid.to_xyz(15, 15, 15);
        REQUIRE_THAT(xyz.x(), Catch::Matchers::WithinAbs(5.0, 1e-6));
        REQUIRE_THAT(xyz.y(), Catch::Matchers::WithinAbs(5.0, 1e-6));
        REQUIRE_THAT(xyz.z(), Catch::Matchers::WithinAbs(5.0, 1e-6));
    }
}

TEST_CASE("Grid::to_x to_y to_z") {
    settings::grid::cell_width = 1.0;
    Limit3D axes(-10, 10, -10, 10, -10, 10);
    Grid grid(axes);

    SECTION("to_x") {
        REQUIRE_THAT(grid.to_x(10), Catch::Matchers::WithinAbs(0.0, 1e-6));
        REQUIRE_THAT(grid.to_x(15), Catch::Matchers::WithinAbs(5.0, 1e-6));
        REQUIRE_THAT(grid.to_x(5), Catch::Matchers::WithinAbs(-5.0, 1e-6));
    }

    SECTION("to_y") {
        REQUIRE_THAT(grid.to_y(10), Catch::Matchers::WithinAbs(0.0, 1e-6));
        REQUIRE_THAT(grid.to_y(15), Catch::Matchers::WithinAbs(5.0, 1e-6));
        REQUIRE_THAT(grid.to_y(5), Catch::Matchers::WithinAbs(-5.0, 1e-6));
    }

    SECTION("to_z") {
        REQUIRE_THAT(grid.to_z(10), Catch::Matchers::WithinAbs(0.0, 1e-6));
        REQUIRE_THAT(grid.to_z(15), Catch::Matchers::WithinAbs(5.0, 1e-6));
        REQUIRE_THAT(grid.to_z(5), Catch::Matchers::WithinAbs(-5.0, 1e-6));
    }
}

TEST_CASE("Grid::roundtrip coordinate conversions") {
    settings::grid::cell_width = 1.0;
    Limit3D axes(-10, 10, -10, 10, -10, 10);
    Grid grid(axes);

    SECTION("to_bins and back to_xyz") {
        Vector3<double> original(3.5, -2.7, 8.1);
        auto bins = grid.to_bins(original);
        auto recovered = grid.to_xyz(bins);
        
        REQUIRE_THAT(recovered.x(), Catch::Matchers::WithinAbs(std::floor(original.x()), 1.0));
        REQUIRE_THAT(recovered.y(), Catch::Matchers::WithinAbs(std::floor(original.y()), 1.0));
        REQUIRE_THAT(recovered.z(), Catch::Matchers::WithinAbs(std::floor(original.z()), 1.0));
    }
}

TEST_CASE("Grid::get_center") {
    settings::grid::cell_width = 1.0;
    Limit3D axes(-10, 10, -10, 10, -10, 10);
    Grid grid(axes);

    auto center = grid.get_center();
    REQUIRE(center.x() == 10);
    REQUIRE(center.y() == 10);
    REQUIRE(center.z() == 10);
}

TEST_CASE("Grid::add Body") {
    settings::grid::cell_width = 1.0;
    Limit3D axes(-10, 10, -10, 10, -10, 10);
    Grid grid(axes);

    SECTION("single atom body") {
        AtomFF atom({0, 0, 0}, form_factor::form_factor_t::C);
        Body body(std::vector<AtomFF>{atom});
        
        REQUIRE(grid.a_members.empty());
        grid.add(body, false);
        REQUIRE(grid.a_members.size() == 1);
        REQUIRE(grid.a_members[0].get_atom() == atom);
    }

    SECTION("multiple atom body") {
        std::vector<AtomFF> atoms = {
            AtomFF({0, 0, 0}, form_factor::form_factor_t::C),
            AtomFF({1, 1, 1}, form_factor::form_factor_t::N),
            AtomFF({-1, -1, -1}, form_factor::form_factor_t::O)
        };
        Body body(atoms);
        
        grid.add(body, false);
        REQUIRE(grid.a_members.size() == 3);
    }
}

TEST_CASE("Grid::add Water") {
    settings::grid::cell_width = 1.0;
    Limit3D axes(-10, 10, -10, 10, -10, 10);
    Grid grid(axes);

    SECTION("single water") {
        Water water({0, 0, 0});
        
        REQUIRE(grid.w_members.empty());
        grid.add(water, false);
        REQUIRE(grid.w_members.size() == 1);
    }

    SECTION("multiple waters") {
        std::vector<Water> waters = {
            Water({0, 0, 0}),
            Water({1, 1, 1}),
            Water({-1, -1, -1})
        };
        
        grid.add(waters, false);
        REQUIRE(grid.w_members.size() == 3);
    }
}

TEST_CASE("Grid::remove Body") {
    settings::grid::cell_width = 1.0;
    Limit3D axes(-10, 10, -10, 10, -10, 10);
    Grid grid(axes);

    std::vector<AtomFF> atoms = {
        AtomFF({0, 0, 0}, form_factor::form_factor_t::C),
        AtomFF({1, 1, 1}, form_factor::form_factor_t::N)
    };
    Body body(atoms);
    
    grid.add(body, false);
    REQUIRE(grid.a_members.size() == 2);
    
    grid.remove(body);
    REQUIRE(grid.a_members.size() == 0);
}

TEST_CASE("Grid::clear_waters") {
    settings::grid::cell_width = 1.0;
    Limit3D axes(-10, 10, -10, 10, -10, 10);
    Grid grid(axes);

    std::vector<Water> waters = {
        Water({0, 0, 0}),
        Water({1, 1, 1}),
        Water({-1, -1, -1})
    };
    
    grid.add(waters, false);
    REQUIRE(grid.w_members.size() == 3);
    
    grid.clear_waters();
    REQUIRE(grid.w_members.empty());
}

TEST_CASE("Grid::get_waters") {
    settings::grid::cell_width = 1.0;
    Limit3D axes(-10, 10, -10, 10, -10, 10);
    Grid grid(axes);

    std::vector<Water> waters = {
        Water({0, 0, 0}),
        Water({1, 1, 1})
    };
    
    grid.add(waters, false);
    auto retrieved_waters = grid.get_waters();
    
    REQUIRE(retrieved_waters.size() == 2);
}

TEST_CASE("Grid::add_volume") {
    settings::grid::cell_width = 1.0;
    Limit3D axes(-10, 10, -10, 10, -10, 10);
    Grid grid(axes);

    int initial_volume = grid.get_volume_bins();
    grid.add_volume(10);
    REQUIRE(grid.get_volume_bins() == initial_volume + 10);
}

TEST_CASE("Grid::index") {
    settings::grid::cell_width = 1.0;
    Limit3D axes(-10, 10, -10, 10, -10, 10);
    Grid grid(axes);

    SECTION("default state is EMPTY") {
        auto state = grid.index(10, 10, 10);
        REQUIRE(state == EMPTY);
    }

    SECTION("after adding atom") {
        AtomFF atom({0, 0, 0}, form_factor::form_factor_t::C);
        Body body(std::vector<AtomFF>{atom});
        grid.add(body, false);
        
        auto bins = grid.to_bins({0, 0, 0});
        auto state = grid.index(bins.x(), bins.y(), bins.z());
        REQUIRE(state == A_CENTER);
    }
}
