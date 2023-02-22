#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rigidbody/RigidBody.h>
#include <rigidbody/transform/RigidTransform.h>
#include <rigidbody/selection/RandomSelect.h>
#include <fitter/HydrationFitter.h>
#include <data/BodySplitter.h>
#include <data/Protein.h>

using namespace rigidbody;

TEST_CASE("can_reuse_fitter", "[files]") {
    setting::general::verbose = false;
    Protein protein_2epe("test/files/2epe.pdb");
    Protein protein_LAR12("data/LAR1-2/LAR1-2.pdb");
    protein_2epe.generate_new_hydration();
    protein_LAR12.generate_new_hydration();

    SECTION("intensity_fitter") {
        fitter::HydrationFitter fitter("test/files/2epe.dat", protein_2epe.get_histogram());
        double chi2 = fitter.fit()->fval;

        fitter.set_scattering_hist(protein_LAR12.get_histogram());
        double _chi2 = fitter.fit()->fval;
        REQUIRE_THAT(chi2, !Catch::Matchers::WithinRel(_chi2));

        fitter.set_scattering_hist(protein_2epe.get_histogram());
        _chi2 = fitter.fit()->fval;
        REQUIRE_THAT(chi2, Catch::Matchers::WithinRel(_chi2));
    }

    SECTION("simple_intensity_fitter") {
        fitter::LinearFitter fitter("test/files/2epe.dat", protein_2epe.get_histogram());
        double chi2 = fitter.fit()->fval;

        fitter.set_scattering_hist(protein_LAR12.get_histogram());
        double _chi2 = fitter.fit()->fval;
        REQUIRE_THAT(chi2, !Catch::Matchers::WithinRel(_chi2));

        fitter.set_scattering_hist(protein_2epe.get_histogram());
        _chi2 = fitter.fit()->fval;
        REQUIRE_THAT(chi2, Catch::Matchers::WithinRel(_chi2));
    }
}

TEST_CASE("can_repeat_fit") {
    setting::general::verbose = false;
    Protein protein("test/files/2epe.pdb");
    fitter::LinearFitter fitter("test/files/2epe.dat", protein.get_histogram());

    protein.generate_new_hydration();
    double chi2 = fitter.fit()->fval;

    for (int i = 0; i < 10; i++) {
        protein.generate_new_hydration();
        double _chi2 = fitter.fit()->fval;
        REQUIRE_THAT(chi2, Catch::Matchers::WithinRel(_chi2));
    }
}

TEST_CASE("rigidbody_opt", "[manual]") {
    setting::general::verbose = false;
    std::vector<int> splits = {9, 99};
    Protein protein(BodySplitter::split("data/LAR1-2/LAR1-2.pdb", splits));
    RigidBody rbody(protein);
    protein.generate_new_hydration();

    fitter::HydrationFitter fitter("data/LAR1-2/LAR1-2.dat", protein.get_histogram());
    fitter.fit()->fval;
}

TEST_CASE("body_selectors") {
    std::vector<Atom> b1 = {Atom(Vector3<double>(-1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, "C", "C", 1)};
    std::vector<Atom> b2 = {Atom(Vector3<double>( 1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>( 1, 1, -1), 1, "C", "C", 1)};
    std::vector<Atom> b3 = {Atom(Vector3<double>(-1, -1,  1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1,  1), 1, "C", "C", 1)};
    std::vector<Atom> b4 = {Atom(Vector3<double>( 1, -1,  1), 1, "C", "C", 1), Atom(Vector3<double>( 1, 1,  1), 1, "C", "C", 1)};
    std::vector<std::vector<Atom>> atoms = {b1, b2, b3, b4};
    Protein protein(atoms, {});

    std::unique_ptr<BodySelectStrategy> strat = std::make_unique<RandomSelect>(protein);
    std::vector<int> count(protein.bodies.size());
    for (unsigned int i = 0; i < 100; i++) {
        unsigned int num = strat->next();
        if (num > count.size()-1) {
            std::cout << "Strategy selected a body outside the allowed range. Number: " << num << std::endl;
            REQUIRE(false);
        } else {
            count[num]++;
        }
    }

    // check that each one was chosen at least 10 times
    REQUIRE(count[0] > 10);
    REQUIRE(count[1] > 10);
    REQUIRE(count[2] > 10);
    REQUIRE(count[3] > 10);
}

TEST_CASE("iteration_step") {
    setting::general::verbose = false;
    auto validate_single_step = [] (Protein& protein) {
        protein.generate_new_hydration();

        // fit the protein
        fitter::HydrationFitter fitter("test/files/2epe.dat", protein.get_histogram());
        auto chi2 = fitter.fit()->fval;

        Body&                 body = protein.bodies.at(0);
        std::shared_ptr<Grid> grid = protein.get_grid();
        Body                  old_body(body);
        Grid                  old_grid = *protein.get_grid();
        Protein               old_protein(protein);

        //####################################//
        //### do one step of rigidbody opt ###//
        //####################################//
        grid->remove(&body);
        body.translate(Vector3<double>(0, 0, 10));              // translate the body
        body.rotate(Vector3<double>(0, 0, 1), 0.1);             // rotate the body
        grid->add(&body);
        protein.generate_new_hydration();                       // generate a new hydration shell
        fitter.set_scattering_hist(protein.get_histogram());    // update the scattering histogram to reflect the new body positions
        auto _chi2 = fitter.fit()->fval;                        // fit the protein
        CHECK_THAT(chi2, !Catch::Matchers::WithinRel(_chi2));   // chi2 should be different

        //######################################//
        //### reset the state of the protein ###//
        //######################################//
        protein.set_grid(old_grid);                             // reset the grid
        body = std::move(old_body);                             // reset the body
        grid = protein.get_grid();                              // reset the grid pointer

        // check that the grid is exactly identical
        auto axes = grid->get_axes();
        auto oldaxes = old_grid.get_axes();
        REQUIRE(axes == oldaxes);
        for (unsigned int i = 0; i < axes.x.bins; i++) {
            for (unsigned int j = 0; j < axes.y.bins; j++) {
                for (unsigned int k = 0; k < axes.z.bins; k++) {
                    CHECK(grid->grid.index(i, j, k) == old_grid.grid.index(i, j, k));
                }
            }
        }
        REQUIRE(*grid == old_grid);

        // check that the atoms are the same
        auto atoms = protein.atoms();
        auto oldatoms = old_protein.atoms();
        REQUIRE(atoms.size() == oldatoms.size());
        for (unsigned int i = 0; i < atoms.size(); i++) {
            CHECK(atoms.at(i).coords == oldatoms.at(i).coords);
        }

        // check that the grid water members are the same
        auto wm = grid->w_members;
        auto oldwm = old_grid.w_members;
        REQUIRE(wm.size() == oldwm.size());

        auto wm_it = wm.begin();
        auto oldwm_it = oldwm.begin();
        while (wm_it != wm.end()) {
            CHECK(*wm_it == *oldwm_it);
            wm_it++;
            oldwm_it++;
        }

        // check that the grid atom members are the same
        auto am = grid->a_members;
        auto oldam = old_grid.a_members;
        REQUIRE(am.size() == oldam.size());

        auto am_it = am.begin();
        auto oldam_it = oldam.begin();
        while (am_it != am.end()) {
            CHECK(*am_it == *oldam_it);
            am_it++;
            oldam_it++;
        }

        // check that there's no water in the grid after a call to clear_hydration
        protein.clear_hydration();
        for (unsigned int i = 0; i < axes.x.bins; i++) {
            for (unsigned int j = 0; j < axes.y.bins; j++) {
                for (unsigned int k = 0; k < axes.z.bins; k++) {
                    REQUIRE(grid->grid.index(i, j, k) != GridObj::H_AREA);
                    REQUIRE(grid->grid.index(i, j, k) != GridObj::H_CENTER);
                }
            }
        }

        //################################################//
        //### check that the protein has been reverted ###//
        //################################################//
        protein.generate_new_hydration();                       // generate a new hydration shell
        fitter.set_scattering_hist(protein.get_histogram());    // update the scattering histogram to reflect the new body positions
        _chi2 = fitter.fit()->fval;                             // fit the protein
        CHECK_THAT(chi2, Catch::Matchers::WithinRel(_chi2));    // chi2 should be the same

        // check that the waters are the same
        auto newwaters = protein.waters();
        auto oldwaters = old_protein.waters();
        REQUIRE(newwaters.size() == oldwaters.size());
        for (unsigned int i = 0; i < newwaters.size(); i++) {
            CHECK(oldwaters.at(i).coords == newwaters.at(i).coords);
        }

        // check that the grid is the same
        for (unsigned int i = 0; i < axes.x.bins; i++) {
            for (unsigned int j = 0; j < axes.y.bins; j++) {
                for (unsigned int k = 0; k < axes.z.bins; k++) {
                    CHECK(grid->grid.index(i, j, k) == old_grid.grid.index(i, j, k));
                }
            }
        }
        REQUIRE(*grid == old_grid);
    };

    SECTION("simple") {
        Atom a1(1, "C", "", "LYS", "", 1, "", Vector3<double>(-1, -1, -1), 1, 0, "C", "");
        Atom a2(2, "C", "", "LYS", "", 1, "", Vector3<double>(-1, -1,  1), 1, 0, "C", "");
        Atom a3(3, "C", "", "LYS", "", 1, "", Vector3<double>(-1,  1, -1), 1, 0, "C", "");
        Atom a4(4, "C", "", "LYS", "", 1, "", Vector3<double>(-1,  1,  1), 1, 0, "C", "");
        Atom a5(5, "C", "", "LYS", "", 1, "", Vector3<double>( 1, -1, -1), 1, 0, "C", "");
        Atom a6(6, "C", "", "LYS", "", 1, "", Vector3<double>( 1, -1,  1), 1, 0, "C", "");
        Atom a7(7, "C", "", "LYS", "", 1, "", Vector3<double>( 1,  1, -1), 1, 0, "C", "");
        Atom a8(8, "C", "", "LYS", "", 1, "", Vector3<double>( 1,  1,  1), 1, 0, "C", "");

        Body b1(std::vector<Atom>{a1, a2});
        Body b2(std::vector<Atom>{a3, a4});
        Body b3(std::vector<Atom>{a5, a6});
        Body b4(std::vector<Atom>{a7, a8});
        std::vector<Body> ap = {b1, b2, b3, b4};
        Protein protein(ap);
        Grid grid(Axis3D(-20, 20, -20, 20, -20, 20, 40), 1);
        grid.add(protein.atoms());
        protein.set_grid(grid);

        validate_single_step(protein);
    }

    SECTION("real data") {
        Protein protein = BodySplitter::split("data/lysozyme/2epe.pdb", {9, 99});
        REQUIRE(protein.bodies.size() == 3);
        validate_single_step(protein);
    }
}

TEST_CASE("iteration_step_old", "[broken]") {
    setting::general::verbose = false;
    Protein protein = BodySplitter::split("data/lysozyme/2epe.pdb", {9, 99});
    REQUIRE(protein.bodies.size() == 3);
    protein.generate_new_hydration();
    std::vector<Water> oldwaters = protein.waters();
    Grid oldgrid = *protein.get_grid();

    // fit the protein
    fitter::HydrationFitter fitter("data/lysozyme/2epe.dat", protein.get_histogram());
    auto chi2 = fitter.fit()->fval;

    // remove the first body
    Body& body = protein.bodies.at(0);
    Body old_body(body);
    std::shared_ptr<Grid> grid = protein.get_grid();

    grid->remove(&body);
    grid->force_expand_volume(); 
    body.translate(Vector3<double>(0, 0, 10));              // translate the body
    body.rotate(Vector3<double>(0, 0, 1), 0.1);             // rotate the body
    grid->add(&body);
    protein.generate_new_hydration();                       // generate a new hydration shell
    fitter.set_scattering_hist(protein.get_histogram());    // update the scattering histogram to reflect the new body positions
    auto _chi2 = fitter.fit()->fval;                        // fit the protein
    CHECK_THAT(chi2, !Catch::Matchers::WithinRel(_chi2));   // chi2 should be different

    // add the body back
    grid->remove(&body);
    grid->force_expand_volume();
    body = old_body;                                        // reset the body
    grid->add(&body);
    protein.generate_new_hydration();                       // generate a new hydration shell
    fitter.set_scattering_hist(protein.get_histogram());    // update the scattering histogram to reflect the new body positions
    _chi2 = fitter.fit()->fval;                             // fit the protein
    CHECK_THAT(chi2, Catch::Matchers::WithinRel(_chi2));    // chi2 should be the same

    // check that the grid is the same
    Grid newgrid = *protein.get_grid();
    REQUIRE(oldgrid == newgrid);

    // check that the waters are the same
    auto newwaters = protein.waters();
    for (unsigned int i = 0; i < newwaters.size(); i++) {
        CHECK(oldwaters.at(i).coords == newwaters.at(i).coords);
    }
}