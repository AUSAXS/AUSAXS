#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hist/DebyeLookupTable.h>
#include <utility/Axis.h>
#include <data/Atom.h>
#include <settings/All.h>

using std::vector;

// Define the hash function for atoms so we can use them in our tests
namespace std {
    template <>
    struct hash<Atom> {
        std::size_t operator()(const Atom& a) const {return a.uid;};
    };
}

TEST_CASE("debye_lookup_table") {
    settings::axes::qmax = 1.001;
    SECTION("default_table") {
        table::DebyeLookupTable::reset();
        double width = settings::axes::distance_bin_width;
        vector<double> d(200/width, 0);
        for (unsigned int i = 1; i < d.size(); i++) {
            d[i] = width*(i+0.5);
        }

        // prepare the q values for the intensity calculations
        Axis debye_axis(settings::axes::qmin, settings::axes::qmax, settings::axes::bins);
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