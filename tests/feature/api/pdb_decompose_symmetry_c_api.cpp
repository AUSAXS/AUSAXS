// Exercises pdb_decompose_symmetry through the C API, as the GUI "view pdb" utility would.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <api/api_pyausaxs.h>

#include <string>

TEST_CASE("pdb_decompose_symmetry: real p2 dimer (SASDJG5)") {
    int status = 0;
    int pdb_id = pdb_read("tests/files/SASDJG5.pdb", &status);
    REQUIRE(status == 0);
    REQUIRE(pdb_id >= 0);

    double *x, *y, *z, rmsd = -1;
    int *copy_index, n_atoms = 0;
    int data_id = pdb_decompose_symmetry(pdb_id, "p2", &x, &y, &z, &copy_index, &n_atoms, &rmsd, &status);
    REQUIRE(status == 0);
    REQUIRE(data_id >= 0);

    // two chains of equal size -> reference (copy 0) + one generated copy (copy 1)
    REQUIRE(n_atoms > 0);
    REQUIRE(n_atoms % 2 == 0);
    int per = n_atoms/2;
    for (int i = 0; i < per; ++i)        {CHECK(copy_index[i] == 0);}
    for (int i = per; i < n_atoms; ++i)  {CHECK(copy_index[i] == 1);}

    // a genuine p2 dimer should decompose with a small residual
    CHECK(0 <= rmsd);
    CHECK(rmsd < 5.0);
}

TEST_CASE("pdb_decompose_symmetry: wrong chain count is reported") {
    int status = 0;
    int pdb_id = pdb_read("tests/files/SASDJG5.pdb", &status);
    REQUIRE(status == 0);

    // SASDJG5 has two chains; c3 would need three
    double *x, *y, *z, rmsd = -1;
    int *copy_index, n_atoms = 0;
    pdb_decompose_symmetry(pdb_id, "c3", &x, &y, &z, &copy_index, &n_atoms, &rmsd, &status);
    CHECK(status != 0);
}
