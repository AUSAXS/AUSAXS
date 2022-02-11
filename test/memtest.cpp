#include "catch2/catch.hpp"

#include "data/Atom.h"
#include "data/Protein.h"
#include "data/BodySplitter.h"

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

TEST_CASE("grid_add", "[memtest]") {
    vector<int> splits = {9, 99};
    Protein protein(BodySplitter::split("data/LAR1-2.pdb", splits));
    RigidBody rbody(protein);
    rbody.protein.generate_new_hydration();
    rbody.protein.clear_grid();
    rbody.protein.generate_new_hydration();
}

TEST_CASE("compact_coordinates", "[memtest]") {
    vector<Atom> a1 = {Atom(1, "C", "", "LYS", "", 1, "", Vector3(-1, -1, -1), 1, 0, "C", "0"),  Atom(2, "C", "", "LYS", "", 1, "", Vector3(-1, 1, -1), 1, 0, "C", "0"),
                       Atom(3, "C", "", "LYS", "", 1, "", Vector3( 1, -1, -1), 1, 0, "C", "0"),  Atom(4, "C", "", "LYS", "", 1, "", Vector3( 1, 1, -1), 1, 0, "C", "0")};
    vector<Atom> a2 = {Atom(5, "C", "", "LYS", "", 1, "", Vector3(-1, -1,  1), 1, 0, "C", "0"),  Atom(6, "C", "", "LYS", "", 1, "", Vector3(-1, 1,  1), 1, 0, "C", "0"),
                       Atom(7, "C", "", "LYS", "", 1, "", Vector3( 1, -1,  1), 1, 0, "C", "0"),  Atom(8, "C", "", "LYS", "", 1, "", Vector3( 1, 1,  1), 1, 0, "C", "0")};

    Body b1(a1);
    Body b2(a2);

    Protein protein({b1, b2});
    RigidBody body(protein);

    vector<double> chi2list = {1100, 900};
    double _chi2 = 1000;
    Parameters params(protein);
    std::shared_ptr<Grid> grid = protein.get_grid();

    for (unsigned int i = 0; i < chi2list.size(); i++) {
        std::cout << "ITERATION " << i << std::endl;
        int body_index = 0;
        Body& body = protein.bodies[body_index];
        Parameter param({1, 1, 1}, 1, 1, 1);

        Body old_body(body);
        Grid old_grid(grid->copy());

        grid->remove(&body);

        Matrix R = Matrix::rotation_matrix(param.alpha, param.beta, param.gamma);
        body.translate(param.dx);
        body.rotate(R);

        grid->add(&body);
        protein.generate_new_hydration();

        double __chi2 = chi2list[i];

        if (__chi2 >= _chi2) {
            protein.bodies[body_index] = old_body;
            protein.generate_new_hydration();
        } else {
            // accept the changes
            _chi2 = __chi2;
            params.update(body.uid, param);
        }
    }
}

TEST_CASE("debug", "[memtest]") {
    vector<Atom> a1 = {Atom(1, "C", "", "LYS", "", 1, "", Vector3(-1, -1, -1), 1, 0, "C", "0"),  Atom(2, "C", "", "LYS", "", 1, "", Vector3(-1, 1, -1), 1, 0, "C", "0"),
                       Atom(3, "C", "", "LYS", "", 1, "", Vector3( 1, -1, -1), 1, 0, "C", "0"),  Atom(4, "C", "", "LYS", "", 1, "", Vector3( 1, 1, -1), 1, 0, "C", "0")};
    vector<Atom> a2 = {Atom(5, "C", "", "LYS", "", 1, "", Vector3(-1, -1,  1), 1, 0, "C", "0"),  Atom(6, "C", "", "LYS", "", 1, "", Vector3(-1, 1,  1), 1, 0, "C", "0"),
                       Atom(7, "C", "", "LYS", "", 1, "", Vector3( 1, -1,  1), 1, 0, "C", "0"),  Atom(8, "C", "", "LYS", "", 1, "", Vector3( 1, 1,  1), 1, 0, "C", "0")};

    // vector<Hetatom> w = {};
    // std::shared_ptr<File> f1 = std::make_shared<File>(a1, w);
    // std::shared_ptr<File> f2 = std::make_shared<File>(a2, w);
    // vector<Atom>& atoms = f1->protein_atoms;

    // {
    //     std::shared_ptr<File> old = std::make_shared<File>(a1, w);
    //     f1 = old;
    // }

    // for (const auto e: f1->protein_atoms) {
    //     std::cout << e.as_pdb() << std::endl;
    // }

    Body b1(a1);
    Body b2(a2);

    {
        Body old_body(b2);
        b1 = old_body;
    }

    for (const auto e : b1.protein_atoms) {
        std::cout << e.as_pdb() << std::endl;
    }
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
