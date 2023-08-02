#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hist/DebyeLookupTable.h>
#include <utility/Axis.h>
#include <data/Atom.h>
#include <settings/All.h>

// Define the hash function for atoms so we can use them in our tests
namespace std {
    template <>
    struct hash<Atom> {
        std::size_t operator()(const Atom& a) const {return a.uid;};
    };
}

TEST_CASE("DebyeLookupTable::DebyeLookupTable") {
    SECTION("vector<double>&, vector<double>&") {
        std::vector<double> d = {1, 2, 3, 4, 5};
        std::vector<double> q = {1, 2, 3, 4, 5};
        table::DebyeLookupTable table(d, q);
        CHECK(table.size_d() == d.size());
        CHECK(table.size_q() == q.size());
    }
}

TEST_CASE("DebyeLookupTable::lookup") {
    std::vector<double> d = {20, 10, 30, 15, 10};
    std::vector<double> q = {0.01, 0.05, 0.1, 0.25, 0.5};
    table::DebyeLookupTable table(d, q);

    auto func = [&] (double q, double d) {return std::sin(q*d)/(q*d);};

    SECTION("double, double") {
        CHECK(table.lookup(q[0], d[0]) == func(q[0], d[0]));
        CHECK(table.lookup(q[1], d[1]) == func(q[1], d[1]));
        CHECK(table.lookup(q[2], d[2]) == func(q[2], d[2]));
        CHECK(table.lookup(q[3], d[3]) == func(q[3], d[3]));
        CHECK(table.lookup(q[4], d[4]) == func(q[4], d[4]));
        CHECK(table.lookup(q[1], d[4]) == func(q[1], d[4]));
        CHECK(table.lookup(q[4], d[2]) == func(q[4], d[2]));
    }

    SECTION("unsigned int, unsigned int") {
        CHECK(table.lookup(1u, 1u) == func(q[1], d[1]));
        CHECK(table.lookup(2u, 2u) == func(q[2], d[2]));
        CHECK(table.lookup(3u, 3u) == func(q[3], d[3]));
        CHECK(table.lookup(4u, 4u) == func(q[4], d[4]));
        CHECK(table.lookup(2u, 5u) == func(q[2], d[5]));
        CHECK(table.lookup(4u, 3u) == func(q[4], d[3]));
    }
}

TEST_CASE("DebyeLookupTable::uses_default_table") {
    SECTION("default table") {
        Axis axis(settings::axes::bins, settings::axes::qmin, settings::axes::qmax);
        std::vector<double> q = axis.as_vector();
        std::vector<double> d = {1, 2, 3, 4, 5};
        table::DebyeLookupTable table(d, q);
        CHECK(table.uses_default_table());
    }
    SECTION("non-default table") {
        std::vector<double> d = {1, 2, 3, 4, 5};
        std::vector<double> q = {1, 2, 3, 4, 5};
        table::DebyeLookupTable table(d, q);
        CHECK(!table.uses_default_table());
    }
}

// size_d & size_q
TEST_CASE("DebyeLookupTable::size") {
    SECTION("(5, 5)") {
        std::vector<double> d = {1, 2, 3, 4, 5};
        std::vector<double> q = {1, 2, 3, 4, 5};
        table::DebyeLookupTable table(d, q);
        CHECK(table.size_d() == d.size());
        CHECK(table.size_q() == q.size());
    }

    SECTION("(5, 7)") {
        std::vector<double> d = {1, 2, 3, 4, 5};
        std::vector<double> q = {1, 2, 3, 4, 5, 6, 7};
        table::DebyeLookupTable table(d, q);
        CHECK(table.size_d() == d.size());
        CHECK(table.size_q() == q.size());
    }

    SECTION("empty") {
        std::vector<double> d = {};
        std::vector<double> q = {};
        table::DebyeLookupTable table(d, q);
        CHECK(table.size_d() == d.size());
        CHECK(table.size_q() == q.size());
    }
}

TEST_CASE("DebyeLookupTable: correct values") {
    settings::axes::qmax = 1.001;
    SECTION("default_table") {
        table::DebyeLookupTable::reset_default_table();
        double width = settings::axes::distance_bin_width;
        std::vector<double> d(200/width, 0);
        for (unsigned int i = 1; i < d.size(); i++) {
            d[i] = width*(i+0.5);
        }

        // prepare the q values for the intensity calculations
        Axis debye_axis(settings::axes::qmin, settings::axes::qmax, settings::axes::bins);
        std::vector<double> q(debye_axis.bins);
        double debye_width = debye_axis.width();
        for (unsigned int i = 0; i < debye_axis.bins; i++) {
            q[i] = debye_axis.min + i*debye_width;
        }

        table::DebyeLookupTable table(q, d);
        REQUIRE(table.uses_default_table());

        CHECK(table.lookup(q[0], d[1]) == sin(q[0]*d[1])/(q[0]*d[1]));
        CHECK(table.lookup(q[1], d[1]) == sin(q[1]*d[1])/(q[1]*d[1]));
        CHECK(table.lookup(q[2], d[1]) == sin(q[2]*d[1])/(q[2]*d[1]));

        CHECK(table.lookup(0u, 0u) == 1);
        CHECK(table.lookup(0u, 1u) == sin(q[0]*d[1])/(q[0]*d[1]));
        CHECK(table.lookup(2u, 1u) == sin(q[2]*d[1])/(q[2]*d[1]));
    }

    SECTION("non_default") {
        std::vector<double> d = {0, 2, 4, 6, 8, 10};
        std::vector<double> q = {0.10, 0.12, 0.14, 0.16, 0.18, 0.20};
        table::DebyeLookupTable table(q, d);
        REQUIRE(!table.uses_default_table());

        CHECK(table.lookup(q[0], d[1]) == sin(q[0]*d[1])/(q[0]*d[1]));
        CHECK(table.lookup(q[1], d[1]) == sin(q[1]*d[1])/(q[1]*d[1]));
        CHECK(table.lookup(q[2], d[1]) == sin(q[2]*d[1])/(q[2]*d[1]));

        CHECK(table.lookup(0u, 0u) == 1);
        CHECK(table.lookup(0u, 1u) == sin(q[0]*d[1])/(q[0]*d[1]));
        CHECK(table.lookup(2u, 1u) == sin(q[2]*d[1])/(q[2]*d[1]));
    }
    settings::axes::qmax = 0.5;
}