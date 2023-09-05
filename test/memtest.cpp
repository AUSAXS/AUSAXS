#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Atom.h>
#include <data/Body.h>
#include <data/Water.h>
#include <data/Protein.h>
#include <data/BodySplitter.h>
#include <math/Matrix.h>
#include <math/MatrixUtils.h>
#include <rigidbody/RigidBody.h>
#include <dataset/SimpleDataset.h>

using std::vector;

TEST_CASE("interpolate_outside_range") {
    std::vector<double> x1, y1;
    for (double xx = 0; xx < 2*M_PI; xx += 0.05) {
        x1.push_back(xx);
        y1.push_back(sin(xx));
    }

    Dataset data1({x1, y1});
    data1.interpolate_y(-1);
    data1.interpolate_y(3*M_PI);
}

/**
 * @brief These tests are meant to be run with valgrind to help identify memory issues. 
 *        I'll expand it as I encounter more leaks/bugs with the software. 
 */
TEST_CASE("atom_assign", "[manual]") {
    vector<Atom> atoms = {Atom(Vector3<double>(-1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, "C", "C", 1),
                          Atom(Vector3<double>( 1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>( 1, 1, -1), 1, "C", "C", 1),
                          Atom(Vector3<double>(-1, -1,  1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1,  1), 1, "C", "C", 1),
                          Atom(Vector3<double>( 1, -1,  1), 1, "C", "C", 1), Atom(Vector3<double>( 1, 1,  1), 1, "C", "C", 1)};

    Body b1(atoms);
    Body b2(atoms);
    b2 = b1;
}

TEST_CASE("body_splitter", "[manual]") {
    vector<int> splits = {9, 99};
    Protein protein = rigidbody::BodySplitter::split("data/LAR1-2/LAR1-2.pdb", splits);
    rigidbody::RigidBody body(protein);
}

TEST_CASE("body_assign", "[manual]") {
    vector<Atom> atoms = {Atom(Vector3<double>(-1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, "C", "C", 1),
                          Atom(Vector3<double>( 1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>( 1, 1, -1), 1, "C", "C", 1),
                          Atom(Vector3<double>(-1, -1,  1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1,  1), 1, "C", "C", 1),
                          Atom(Vector3<double>( 1, -1,  1), 1, "C", "C", 1), Atom(Vector3<double>( 1, 1,  1), 1, "C", "C", 1)};

    Body b1(atoms);
    Body b2(atoms);

    b2 = b1;
}

#include <plots/PlotDataset.h>
TEST_CASE("memtest_rebin", "[files]") {
    SimpleDataset data("data/SHOC2/SHOC2.dat");
    SimpleDataset data_unbinned = data;
    // data.rebin();
    // data.save("temp/dataset/rebin.dat");

    // plots::PlotDataset plot(data_unbinned);
    // plot.plot(data);
    // plot.save("temp/dataset_rebin.pdf");
}

TEST_CASE("vector_assign", "[manual]") {
    vector<int> v = {1, 2, 3, 4};
    int& v1 = v[0];
    v.assign({0, 3});
    std::cout << v1 << std::endl;
}

TEST_CASE("grid_add", "[manual]") {
    vector<int> splits = {9, 99};
    rigidbody::RigidBody rbody = rigidbody::BodySplitter::split("data/LAR1-2/LAR1-2.pdb", splits);
    rbody.generate_new_hydration();
    rbody.clear_grid();
    rbody.generate_new_hydration();
}

TEST_CASE("debug", "[manual]") {
    vector<Atom> a1 = {Atom(1, "C", "", "LYS", "", 1, "", Vector3<double>(-1, -1, -1), 1, 0, "C", "0"),  Atom(2, "C", "", "LYS", "", 1, "", Vector3<double>(-1, 1, -1), 1, 0, "C", "0"),
                       Atom(3, "C", "", "LYS", "", 1, "", Vector3<double>( 1, -1, -1), 1, 0, "C", "0"),  Atom(4, "C", "", "LYS", "", 1, "", Vector3<double>( 1, 1, -1), 1, 0, "C", "0")};
    vector<Atom> a2 = {Atom(5, "C", "", "LYS", "", 1, "", Vector3<double>(-1, -1,  1), 1, 0, "C", "0"),  Atom(6, "C", "", "LYS", "", 1, "", Vector3<double>(-1, 1,  1), 1, 0, "C", "0"),
                       Atom(7, "C", "", "LYS", "", 1, "", Vector3<double>( 1, -1,  1), 1, 0, "C", "0"),  Atom(8, "C", "", "LYS", "", 1, "", Vector3<double>( 1, 1,  1), 1, 0, "C", "0")};

    Body b1(a1);
    Body b2(a2);

    {
        Body old_body(b2);
        b1 = old_body;
    }

    for (const auto& e : b1.get_atoms()) {
        std::cout << e.as_pdb() << std::endl;
    }
}

TEST_CASE("debug2", "[manual]") {
    vector<Atom> a1 = {Atom(1, "C", "", "LYS", "", 1, "", Vector3<double>(-1, -1, -1), 1, 0, "C", "0"),  Atom(2, "C", "", "LYS", "", 1, "", Vector3<double>(-1, 1, -1), 1, 0, "C", "0"),
                       Atom(3, "C", "", "LYS", "", 1, "", Vector3<double>( 1, -1, -1), 1, 0, "C", "0"),  Atom(4, "C", "", "LYS", "", 1, "", Vector3<double>( 1, 1, -1), 1, 0, "C", "0")};
    vector<Atom> a2 = {Atom(5, "C", "", "LYS", "", 1, "", Vector3<double>(-1, -1,  1), 1, 0, "C", "0"),  Atom(6, "C", "", "LYS", "", 1, "", Vector3<double>(-1, 1,  1), 1, 0, "C", "0"),
                       Atom(7, "C", "", "LYS", "", 1, "", Vector3<double>( 1, -1,  1), 1, 0, "C", "0"),  Atom(8, "C", "", "LYS", "", 1, "", Vector3<double>( 1, 1,  1), 1, 0, "C", "0")};

    vector<Water> w = {};
    std::shared_ptr<ProteinFile> f1 = std::make_shared<ProteinFile>(a1, w);
    std::shared_ptr<ProteinFile> f2 = std::make_shared<ProteinFile>(a2, w);
    // vector<Atom>& atoms = f1->protein_atoms;

    {
        std::shared_ptr<ProteinFile> old = std::make_shared<ProteinFile>(a1, w);
        f1 = old;
    }

    for (const auto& e: f1->protein_atoms) {
        std::cout << e.as_pdb() << std::endl;
    }
}

TEST_CASE("file_assign", "[manual]") {
    vector<Atom> atoms = {Atom(Vector3<double>(-1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, "C", "C", 1),
                          Atom(Vector3<double>( 1, -1, -1), 1, "C", "C", 1), Atom(Vector3<double>( 1, 1, -1), 1, "C", "C", 1),
                          Atom(Vector3<double>(-1, -1,  1), 1, "C", "C", 1), Atom(Vector3<double>(-1, 1,  1), 1, "C", "C", 1),
                          Atom(Vector3<double>( 1, -1,  1), 1, "C", "C", 1), Atom(Vector3<double>( 1, 1,  1), 1, "C", "C", 1)};

    vector<Water> waters;

    std::shared_ptr<ProteinFile> file2 = std::make_shared<ProteinFile>();
    std::shared_ptr<ProteinFile> file1 = std::make_shared<ProteinFile>(atoms, waters);

    // MEMORY CORRUPT
    // vector<Atom>& file2a = file2->protein_atoms;
    // file2 = std::make_shared<File>(*file1);
    // file2a = file2->protein_atoms;
}
