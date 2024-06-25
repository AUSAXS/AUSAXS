#include <catch2/catch_test_macros.hpp>

#include <grid/Grid.h>
#include <grid/detail/GridSurfaceDetection.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <data/record/Atom.h>
#include <settings/All.h>
 
#include <vector>
#include <string>

using namespace data;
using namespace data::record;

class GridDebug : public grid::Grid {
    public: 
        using Grid::Grid;

		double get_atomic_radius(constants::atom_t) const override {return ra;}
		double get_hydration_radius() const override {return rh;}
        void set_atomic_radius(double ra) {this->ra = ra;}
        void set_hydration_radius(double rh) {this->rh = rh;}

        static void generate_debug_grid(Molecule& protein) {
            settings::grid::min_bins = 50;
            auto grid = std::make_unique<GridDebug>(protein.get_bodies());
            grid->set_atomic_radius(0);
            protein.set_grid(std::move(grid));
        }
    
    private:
        double ra = 0, rh = 0;
};

// Check that the GridSurfaceDetection::detect_atoms function works as expected.
TEST_CASE("GridSurfaceDetection::detect_atoms") {
    settings::molecule::use_effective_charge = false;
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = true;

    SECTION("Single with radius") {
        // test a single atom with radius sqrt(3), such that it will occupy a 3x3x3 cube, with everything but the center being 'surface'. 
        settings::grid::rvol = std::sqrt(3)+1e-3;

        Molecule protein({Atom({0, 0, 0}, 1, constants::atom_t::C, "C", 1)});
        GridDebug::generate_debug_grid(protein);
        auto vol = protein.get_grid()->generate_excluded_volume(true);

        CHECK(vol.interior.size() == 1);
        CHECK(vol.surface.size() == 26);
    }

    SECTION("Two atoms with radius") {
        // test a simple 2 atom system with 1Å radii. Now the interior 2 voxels should be 'interior', with the rest 'surface'.
        settings::grid::rvol = std::sqrt(2)+1e-3;
        std::vector<Atom> atoms = {
            Atom({0, 0, 0}, 1, constants::atom_t::C, "C", 1),
            Atom({1, 0, 0}, 1, constants::atom_t::C, "C", 1)
        };

        Molecule protein(atoms);
        GridDebug::generate_debug_grid(protein);
        auto vol = protein.get_grid()->generate_excluded_volume(true);

        CHECK(vol.interior.size() == 2);
        CHECK(vol.surface.size() == 26);
    }

    SECTION("2x2x2 interior simple") {
        // test a simple 2x2x2 system with 2Å radii. 
        settings::grid::rvol = 2;
        std::vector<Atom> atoms;
        for (double x = 0; x <= 1; x+=1) {
            for (double y = 0; y <= 1; y+=1) {
                for (double z = 0; z <= 1; z+=1) {
                    atoms.push_back(Atom({x, y, z}, 1, constants::atom_t::C, "C", 1));
                }
            }
        }
        REQUIRE(atoms.size() == 8);

        Molecule protein(atoms);
        GridDebug::generate_debug_grid(protein);
        auto vol = protein.get_grid()->generate_excluded_volume(true);

        CHECK(vol.interior.size() == 32);
        CHECK(vol.surface.size()  == 56);
    }

    SECTION("2x2x2 interior advanced") {
        // test a simple 2x2x2 system with 1Å radii. Now the interior 2x2x2 cube should be 'interior', with the rest 'surface'.
        settings::grid::rvol = std::sqrt(2)+1e-3;
        std::vector<Atom> atoms;
        for (double x = 0; x <= 1; x+=1) {
            for (double y = 0; y <= 1; y+=1) {
                for (double z = 0; z <= 1; z+=1) {
                    atoms.push_back(Atom({x, y, z}, 1, constants::atom_t::C, "C", 1));
                }
            }
        }
        REQUIRE(atoms.size() == 8);

        Molecule protein(atoms);
        GridDebug::generate_debug_grid(protein);
        auto vol = protein.get_grid()->generate_excluded_volume(true);

        CHECK(vol.interior.size() == 8);
        CHECK(vol.surface.size() == 48);
    }

    SECTION("3x3x3") {
        // test a simple 3x3x3 system with 1Å radii. Now the interior 3x3x3 cube should be 'interior', with the rest 'surface'.
        settings::grid::rvol = 2;
        std::vector<Atom> atoms;
        for (double x = -1; x <= 1; x+=1) {
            for (double y = -1; y <= 1; y+=1) {
                for (double z = -1; z <= 1; z+=1) {
                    atoms.push_back(Atom({x, y, z}, 1, constants::atom_t::C, "C", 1));
                }
            }
        }
        REQUIRE(atoms.size() == 27);

        Molecule protein(atoms);
        GridDebug::generate_debug_grid(protein);
        auto vol = protein.get_grid()->generate_excluded_volume(true);

        CHECK(vol.interior.size() == 81);
        CHECK( vol.surface.size() == 98);
    }

    SECTION("larger radius, single") {
        settings::grid::rvol = 2;
        std::vector<Atom> atoms = {Atom({0, 0, 0}, 1, constants::atom_t::C, "C", 1)};

        Molecule protein(atoms);
        GridDebug::generate_debug_grid(protein);
        auto vol = protein.get_grid()->generate_excluded_volume(true);

        CHECK(vol.interior.size() == 1);
        CHECK(vol.surface.size() == 32);
    }

    SECTION("much larger radius, single") {
        settings::grid::rvol = 3;
        std::vector<Atom> atoms = {Atom({0, 0, 0}, 1, constants::atom_t::C, "C", 1)};

        Molecule protein(atoms);
        GridDebug::generate_debug_grid(protein);
        auto vol = protein.get_grid()->generate_excluded_volume(true);

        CHECK(vol.interior.size() == 27);
        CHECK(vol.surface.size()  == 96);
    }

    SECTION("larger radius, cube") {
        settings::grid::rvol = 2;
        std::vector<Atom> atoms;
        for (double x = -1; x <= 1; x+=1) {
            for (double y = -1; y <= 1; y+=1) {
                for (double z = -1; z <= 1; z+=1) {
                    atoms.push_back(Atom({x, y, z}, 1, constants::atom_t::C, "C", 1));
                }
            }
        }
        REQUIRE(atoms.size() == 27);

        Molecule protein(atoms);
        GridDebug::generate_debug_grid(protein);
        auto vol = protein.get_grid()->generate_excluded_volume(true);

        CHECK(vol.interior.size() == 81);
        CHECK(vol.surface.size() == 98);
    }

    SECTION("larger radius, larger cube") {
        settings::grid::rvol = 2;
        std::vector<Atom> atoms = {
            Atom({-2, 0, 0}, 1, constants::atom_t::C, "C", 1),
            Atom({2, 0, 0}, 1, constants::atom_t::C, "C", 1),
            Atom({0, -2, 0}, 1, constants::atom_t::C, "C", 1),
            Atom({0, 2, 0}, 1, constants::atom_t::C, "C", 1)
        };

        Molecule protein(atoms);
        GridDebug::generate_debug_grid(protein);
        auto vol = protein.get_grid()->generate_excluded_volume(true);

        CHECK(vol.interior.size() == 21);
        CHECK(vol.surface.size() == 92);
    }

    // SECTION("6lyz_exv") {
    //     settings::grid::save_exv = true;
    //     settings::general::output = "temp/tests/grid/";
    //     settings::grid::rvol = 3;
    //     Molecule protein("tests/files/6lyz_exv.pdb");
    //     // protein.get_body(0).get_atoms().resize(20);
    //     GridDebug::generate_debug_grid(protein);
    //     auto vol = protein.get_grid()->generate_excluded_volume(true);
    // }
}

TEST_CASE("GridSurfaceDetection::detect_voxels") {
    settings::molecule::use_effective_charge = false;
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = true;

    settings::grid::save_exv = true;
    settings::general::output = "temp/tests/grid/";

    settings::grid::rvol = 2.5;
    std::vector<Atom> atoms = {
        Atom({-4,  0,  0}, 1, constants::atom_t::dummy, "C", 1),
        Atom({ 4,  0,  0}, 1, constants::atom_t::dummy, "C", 1),
        Atom({ 0, -4,  0}, 1, constants::atom_t::dummy, "C", 1),
        Atom({ 0,  4,  0}, 1, constants::atom_t::dummy, "C", 1),
        Atom({ 0,  0,  4}, 1, constants::atom_t::dummy, "C", 1),
        Atom({ 0,  0, -4}, 1, constants::atom_t::dummy, "C", 1)
    };

    Molecule protein(atoms);
    GridDebug::generate_debug_grid(protein);
    static_cast<GridDebug*>(protein.get_grid())->set_hydration_radius(3);
    protein.generate_new_hydration();
    auto vol = protein.get_grid()->generate_excluded_volume(true);
    protein.get_grid()->save("temp/tests/grid/exv_full.pdb");

    CHECK(vol.interior.size() == 21);
    CHECK(vol.surface.size() == 92);
}