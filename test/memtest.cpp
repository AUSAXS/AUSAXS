#include <catch2/catch_all.hpp>

#include <data/Atom.h>
#include <data/Protein.h>
#include <data/BodySplitter.h>
#include <math/Matrix.h>

/**
 * @brief These tests are meant to be run with valgrind to help identify memory issues. 
 *        I'll expand it as I encounter more leaks/bugs with the tool. 
 */
TEST_CASE("atom_assign", "[memtest]") {
    vector<Atom> atoms = {Atom(Vector3(-1, -1, -1), 1, "C", "C", 1), Atom(Vector3(-1, 1, -1), 1, "C", "C", 1),
                          Atom(Vector3( 1, -1, -1), 1, "C", "C", 1), Atom(Vector3( 1, 1, -1), 1, "C", "C", 1),
                          Atom(Vector3(-1, -1,  1), 1, "C", "C", 1), Atom(Vector3(-1, 1,  1), 1, "C", "C", 1),
                          Atom(Vector3( 1, -1,  1), 1, "C", "C", 1), Atom(Vector3( 1, 1,  1), 1, "C", "C", 1)};

    Body b1(atoms);
    Body b2(atoms);
    b2 = b1;
}

TEST_CASE("body_splitter", "[memtest]") {
    vector<int> splits = {9, 99};
    Protein protein = BodySplitter::split("data/LAR1-2.pdb", splits);
    RigidBody body(protein);
}

TEST_CASE("body_assign", "[memtest]") {
    vector<Atom> atoms = {Atom(Vector3(-1, -1, -1), 1, "C", "C", 1), Atom(Vector3(-1, 1, -1), 1, "C", "C", 1),
                          Atom(Vector3( 1, -1, -1), 1, "C", "C", 1), Atom(Vector3( 1, 1, -1), 1, "C", "C", 1),
                          Atom(Vector3(-1, -1,  1), 1, "C", "C", 1), Atom(Vector3(-1, 1,  1), 1, "C", "C", 1),
                          Atom(Vector3( 1, -1,  1), 1, "C", "C", 1), Atom(Vector3( 1, 1,  1), 1, "C", "C", 1)};

    Body b1(atoms);
    Body b2(atoms);

    b2 = b1;
}

TEST_CASE("vector_assign", "[memtest]") {
    vector<int> v = {1, 2, 3, 4};
    int& v1 = v[0];
    v.assign({0, 3});
    std::cout << v1 << std::endl;
}

TEST_CASE("grid_add", "[memtest]") {
    vector<int> splits = {9, 99};
    Protein protein(BodySplitter::split("data/LAR1-2.pdb", splits));
    RigidBody rbody(protein);
    rbody.protein.generate_new_hydration();
    rbody.protein.clear_grid();
    rbody.protein.generate_new_hydration();
}

#include "rigidbody/RandomSelect.h"
#include "rigidbody/SimpleParameterGeneration.h"
TEST_CASE("compact_coordinates", "[memtest]") {
    vector<Atom> a1 = {Atom(1, "C", "", "LYS", "", 1, "", Vector3(-1, -1, -1), 1, 0, "C", "0"),  Atom(2, "C", "", "LYS", "", 1, "", Vector3(-1, 1, -1), 1, 0, "C", "0"),
                       Atom(3, "C", "", "LYS", "", 1, "", Vector3( 1, -1, -1), 1, 0, "C", "0"),  Atom(4, "C", "", "LYS", "", 1, "", Vector3( 1, 1, -1), 1, 0, "C", "0")};
    vector<Atom> a2 = {Atom(5, "C", "", "LYS", "", 1, "", Vector3(-1, -1,  1), 1, 0, "C", "0"),  Atom(6, "C", "", "LYS", "", 1, "", Vector3(-1, 1,  1), 1, 0, "C", "0"),
                       Atom(7, "C", "", "LYS", "", 1, "", Vector3( 1, -1,  1), 1, 0, "C", "0"),  Atom(8, "C", "", "LYS", "", 1, "", Vector3( 1, 1,  1), 1, 0, "C", "0")};

    Body b1(a1);
    Body b2(a2);

    Protein protein = BodySplitter::split("data/LAR1-2.pdb", {9, 99});
    RigidBody rigidbody(protein);

    Parameters params(protein);
    // std::shared_ptr<Grid> grid = protein.get_grid();

    rigidbody.generate_new_hydration();
    SimpleIntensityFitter fitter("data/LAR1-2.RSR", protein.get_histogram());
    double _chi2 = fitter.fit()->chi2;
    std::cout << "Initial chi2: " << _chi2 << std::endl;

    RandomSelect bodyselector(rigidbody.protein);
    SimpleParameterGeneration parameter_generator(100, 5, 0.3);

    for (unsigned int i = 0; i < 100; i++) {
        // select a body to be modified this iteration
        int body_index = bodyselector.next();
        Body& body = rigidbody.protein.bodies[body_index];
        Parameter param = parameter_generator.next();

        Body old_body(body);
        // grid->remove(&body);

        // update the body to reflect the new params
        Matrix R = Matrix<double>::rotation_matrix(param.alpha, param.beta, param.gamma);
        body.translate(param.dx);
        body.rotate(R);
 
        // add the body to the grid again
        // grid->add(&body);
        std::shared_ptr<Grid> grid = protein.create_grid();
        rigidbody.generate_new_hydration();

        // calculate the new chi2
        auto h = rigidbody.protein.get_histogram();
        fitter.set_scattering_hist(h);
        double __chi2 = fitter.fit()->chi2;

        std::cout << "new chi2: " << __chi2 << std::endl;
        // if the old configuration was better
        if (__chi2 >= _chi2) {
            std::cout << "REROLLING CHANGES" << std::endl;
            body = old_body;
        } else {
            std::cout << "KEEPING CHANGES" << std::endl;
            // accept the changes
            _chi2 = __chi2;
            params.update(body.uid, param);
        }

        grid = protein.create_grid();
        rigidbody.generate_new_hydration();
        fitter.set_scattering_hist(protein.get_histogram());
        double ___chi2 = fitter.fit()->chi2;
        std::cout << "\tchi2 is now " << ___chi2 << std::endl;
    }
}

TEST_CASE("debug", "[memtest]") {
    vector<Atom> a1 = {Atom(1, "C", "", "LYS", "", 1, "", Vector3(-1, -1, -1), 1, 0, "C", "0"),  Atom(2, "C", "", "LYS", "", 1, "", Vector3(-1, 1, -1), 1, 0, "C", "0"),
                       Atom(3, "C", "", "LYS", "", 1, "", Vector3( 1, -1, -1), 1, 0, "C", "0"),  Atom(4, "C", "", "LYS", "", 1, "", Vector3( 1, 1, -1), 1, 0, "C", "0")};
    vector<Atom> a2 = {Atom(5, "C", "", "LYS", "", 1, "", Vector3(-1, -1,  1), 1, 0, "C", "0"),  Atom(6, "C", "", "LYS", "", 1, "", Vector3(-1, 1,  1), 1, 0, "C", "0"),
                       Atom(7, "C", "", "LYS", "", 1, "", Vector3( 1, -1,  1), 1, 0, "C", "0"),  Atom(8, "C", "", "LYS", "", 1, "", Vector3( 1, 1,  1), 1, 0, "C", "0")};

    Body b1(a1);
    Body b2(a2);

    {
        Body old_body(b2);
        b1 = old_body;
    }

    for (const auto& e : b1.protein_atoms) {
        std::cout << e.as_pdb() << std::endl;
    }
}

TEST_CASE("debug2", "[memtest]") {
    vector<Atom> a1 = {Atom(1, "C", "", "LYS", "", 1, "", Vector3(-1, -1, -1), 1, 0, "C", "0"),  Atom(2, "C", "", "LYS", "", 1, "", Vector3(-1, 1, -1), 1, 0, "C", "0"),
                       Atom(3, "C", "", "LYS", "", 1, "", Vector3( 1, -1, -1), 1, 0, "C", "0"),  Atom(4, "C", "", "LYS", "", 1, "", Vector3( 1, 1, -1), 1, 0, "C", "0")};
    vector<Atom> a2 = {Atom(5, "C", "", "LYS", "", 1, "", Vector3(-1, -1,  1), 1, 0, "C", "0"),  Atom(6, "C", "", "LYS", "", 1, "", Vector3(-1, 1,  1), 1, 0, "C", "0"),
                       Atom(7, "C", "", "LYS", "", 1, "", Vector3( 1, -1,  1), 1, 0, "C", "0"),  Atom(8, "C", "", "LYS", "", 1, "", Vector3( 1, 1,  1), 1, 0, "C", "0")};

    vector<Hetatom> w = {};
    std::shared_ptr<File> f1 = std::make_shared<File>(a1, w);
    std::shared_ptr<File> f2 = std::make_shared<File>(a2, w);
    // vector<Atom>& atoms = f1->protein_atoms;

    {
        std::shared_ptr<File> old = std::make_shared<File>(a1, w);
        f1 = old;
    }

    for (const auto& e: f1->protein_atoms) {
        std::cout << e.as_pdb() << std::endl;
    }
}

TEST_CASE("test", "[memtest]") {
    vector<double> v1(0);
    vector<double>& vref = v1; 
    vector<double> v2(vref);

    std::cout << v2.size() << std::endl;
}

TEST_CASE("file_assign", "[memtest]") {
    vector<Atom> atoms = {Atom(Vector3(-1, -1, -1), 1, "C", "C", 1), Atom(Vector3(-1, 1, -1), 1, "C", "C", 1),
                          Atom(Vector3( 1, -1, -1), 1, "C", "C", 1), Atom(Vector3( 1, 1, -1), 1, "C", "C", 1),
                          Atom(Vector3(-1, -1,  1), 1, "C", "C", 1), Atom(Vector3(-1, 1,  1), 1, "C", "C", 1),
                          Atom(Vector3( 1, -1,  1), 1, "C", "C", 1), Atom(Vector3( 1, 1,  1), 1, "C", "C", 1)};

    vector<Hetatom> waters;

    std::shared_ptr<File> file2 = std::make_shared<File>();
    std::shared_ptr<File> file1 = std::make_shared<File>(atoms, waters);

    // MEMORY CORRUPT
    // vector<Atom>& file2a = file2->protein_atoms;
    // file2 = std::make_shared<File>(*file1);
    // file2a = file2->protein_atoms;
}
