#include <Table.h>
#include <data/Atom.h>

#include "catch2/catch.hpp"

TEST_CASE("lookup tables", "[table]") {
    SECTION("integer table") {
        vector<int> rows = {1, 2, 3};
        vector<int> cols = {1, 2, 3};
        LookupTable<int, int> table(rows, cols);

        for (int i = 0; i < rows.size(); i++) {
            for (int j = 0; j < cols.size(); j++) {
                table.assign(i, j, i*j);
            }
        }

        CHECK(table.lookup(1, 1) == 1);
        CHECK(table.lookup(1, 2) == 2);
        CHECK(table.lookup(1, 3) == 3);
        CHECK(table.lookup(2, 1) == 2);
        CHECK(table.lookup(2, 2) == 4);
        CHECK(table.lookup(2, 3) == 6);
        CHECK(table.lookup(3, 1) == 3);
        CHECK(table.lookup(3, 2) == 6);
        CHECK(table.lookup(3, 3) == 9);
    }

    SECTION("double table") {
        vector<double> rows = {0.35, -2.5, 10.2};
        vector<double> cols = {-3.1, -5.05, 6.7};
        LookupTable<double, double> table(rows, cols);

        for (int i = 0; i < rows.size(); i++) {
            for (int j = 0; j < cols.size(); j++) {
                table.assign(rows[i], cols[j], rows[i]+cols[j]);
            }
        }

        CHECK(table.lookup(1, 1) == 0.35 - 3.1);
        CHECK(table.lookup(1, 2) == 0.35 - 5.05);
        CHECK(table.lookup(1, 3) == 0.35 + 6.7);
        CHECK(table.lookup(2, 1) == -2.5 - 3.1);
        CHECK(table.lookup(2, 2) == -2.5 - 5.05);
        CHECK(table.lookup(2, 3) == -2.5 + 6.7);
        CHECK(table.lookup(3, 1) == 10.2 - 3.1);
        CHECK(table.lookup(3, 2) == 10.2 - 5.05);
        CHECK(table.lookup(3, 3) == 10.2 + 6.7);        
    }

    SECTION("object table") {
        vector<Atom> atoms1 = {Atom(Vector3(-1, -1, -1), 1, "C", "C", 1), Atom(Vector3(-1, 1, -1), 1, "C", "C", 1), Atom(Vector3(1, -1, -1), 1, "C", "C", 1)};
        vector<Atom> atoms2 = {Atom(Vector3(-1, -1, 1), 1, "C", "C", 1), Atom(Vector3(-1, 1, 1), 1, "C", "C", 1)};

        LookupTable<Atom, Atom> table(atoms1, atoms2);

        for (int i = 0; i < atoms1.size(); i++) {
            for (int j = 0; j < atoms2.size(); j++) {
                table.assign(atoms1[i], atoms2[j], i+j);
            }
        }

        CHECK(table.lookup(atoms1[0], atoms2[0]) == 0);
        CHECK(table.lookup(atoms1[0], atoms2[1]) == 1);
        CHECK(table.lookup(atoms1[1], atoms2[0]) == 1);
        CHECK(table.lookup(atoms1[1], atoms2[1]) == 2);
        CHECK(table.lookup(atoms1[2], atoms2[0]) == 2);
        CHECK(table.lookup(atoms1[2], atoms2[1]) == 3);
    }
}

TEST_CASE("debye_lookup_table", "[table]") {
    
}