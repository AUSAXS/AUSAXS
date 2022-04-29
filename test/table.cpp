#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <histogram/DebyeLookupTable.h>
#include <data/Atom.h>

// Define the hash function for atoms so we can use them in our tests
namespace std {
    template <>
    struct hash<Atom> {
        std::size_t operator()(const Atom& a) const {return a.uid;};
    };
}

TEST_CASE("lookup_tables", "[table]") {
    SECTION("integer table") {
        vector<unsigned int> rows = {1, 2, 3};
        vector<unsigned int> cols = {1, 2, 3};
        LookupTable<unsigned int, unsigned int> table(rows, cols);

        for (unsigned int r : rows) {
            for (unsigned int c : cols) {
                table.assign(r, c, r*c);
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

        for (double r : rows) {
            for (double c : cols) {
                table.assign(r, c, r+c);
            }
        }

        CHECK(table.lookup(rows[0], cols[0]) == rows[0] + cols[0]);
        CHECK(table.lookup(rows[0], cols[1]) == rows[0] + cols[1]);
        CHECK(table.lookup(rows[0], cols[2]) == rows[0] + cols[2]);
        CHECK(table.lookup(rows[1], cols[0]) == rows[1] + cols[0]);
        CHECK(table.lookup(rows[1], cols[1]) == rows[1] + cols[1]);
        CHECK(table.lookup(rows[1], cols[2]) == rows[1] + cols[2]);
        CHECK(table.lookup(rows[2], cols[0]) == rows[2] + cols[0]);
        CHECK(table.lookup(rows[2], cols[1]) == rows[2] + cols[1]);
        CHECK(table.lookup(rows[2], cols[2]) == rows[2] + cols[2]);        
    }

    SECTION("object table") {
        vector<Atom> atoms1 = {Atom(Vector3(-1, -1, -1), 1, "C", "C", 1), Atom(Vector3(-1, 1, -1), 1, "C", "C", 1), Atom(Vector3(1, -1, -1), 1, "C", "C", 1)};
        vector<Atom> atoms2 = {Atom(Vector3(-1, -1, 1), 1, "C", "C", 1), Atom(Vector3(-1, 1, 1), 1, "C", "C", 1)};

        LookupTable<Atom, Atom> table(atoms1, atoms2);

        for (unsigned int i = 0; i < atoms1.size(); i++) {
            for (unsigned int j = 0; j < atoms2.size(); j++) {
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
    SECTION("default_table") {
        double width = setting::axes::scattering_intensity_plot_binned_width;
        vector<double> d(200/width, 0);
        for (unsigned int i = 1; i < d.size(); i++) {
            d[i] = width*(i+0.5);
        }

        // prepare the q values for the intensity calculations
        const Axis& debye_axis = setting::axes::scattering_intensity_plot_axis;
        vector<double> q(debye_axis.bins);
        double debye_width = debye_axis.width();
        for (unsigned int i = 0; i < debye_axis.bins; i++) {
            q[i] = debye_axis.min + i*debye_width;
        }

        table::DebyeLookupTable table(q, d);
        REQUIRE(table.uses_default_table());

        CHECK(table.lookup(q[0], d[1]) == sin(q[0]*d[1])/(q[0]*d[1]));
        CHECK(table.lookup(q[1], d[1]) == sin(q[1]*d[1])/(q[1]*d[1]));
        CHECK(table.lookup(q[2], d[1]) == sin(q[2]*d[1])/(q[2]*d[1]));

        CHECK(table.lookup(0, 0) == 1);
        CHECK(table.lookup(0, 1) == sin(q[0]*d[1])/(q[0]*d[1]));
        CHECK(table.lookup(2, 1) == sin(q[2]*d[1])/(q[2]*d[1]));
    }

    SECTION("non_default") {
        vector<double> d = {0, 2, 4, 6, 8, 10};
        vector<double> q = {0.10, 0.12, 0.14, 0.16, 0.18, 0.20};
        table::DebyeLookupTable table(q, d);
        REQUIRE(!table.uses_default_table());

        CHECK(table.lookup(q[0], d[1]) == sin(q[0]*d[1])/(q[0]*d[1]));
        CHECK(table.lookup(q[1], d[1]) == sin(q[1]*d[1])/(q[1]*d[1]));
        CHECK(table.lookup(q[2], d[1]) == sin(q[2]*d[1])/(q[2]*d[1]));

        CHECK(table.lookup(0, 0) == 1);
        CHECK(table.lookup(0, 1) == sin(q[0]*d[1])/(q[0]*d[1]));
        CHECK(table.lookup(2, 1) == sin(q[2]*d[1])/(q[2]*d[1]));
    }
}