#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <table/VectorDebyeTable.h>

#include <cmath>

using namespace ausaxs;

static auto& q = constants::axes::q_vals;
static auto& d = constants::axes::d_vals;

TEST_CASE("VectorDebyeTable: correct zero values") {
    table::VectorDebyeTable debye_table(std::vector<double>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, constants::axes::q_vals);
    SECTION("d = 0") {
        for (int i = 0; i < static_cast<int>(q.size()); ++i) {
            CHECK(debye_table.lookup(i, 0) == 1);
        }
    }
}

TEST_CASE("VectorDebyeTable: correct values") {
    table::VectorDebyeTable debye_table(constants::axes::d_vals, constants::axes::q_vals);
    SECTION("default_table") {
        REQUIRE(debye_table.size_d() == d.size());
        REQUIRE(debye_table.size_q() == q.size());
        for (int i = 0; i < static_cast<int>(q.size()); ++i) {
            for (int j = 0; j < static_cast<int>(d.size()); ++j) {
                double qd = q[i]*d[j];
                if (qd < 1e-3) {
                    CHECK_THAT(debye_table.lookup(i, j), Catch::Matchers::WithinAbs(1 - qd*qd/6 + qd*qd*qd*qd/120, 1e-3)); 
                    continue;
                }
                CHECK_THAT(debye_table.lookup(i, j), Catch::Matchers::WithinAbs(std::sin(qd)/qd, 1e-3));
            }
        }
    }
}

TEST_CASE("VectorDebyeTable: iterators") {
    table::VectorDebyeTable debye_table(std::vector<double>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, constants::axes::q_vals);
    for (int i = 0; i < 10; ++i) {
        int j = 0;
        for (auto* it = debye_table.begin(i); it != debye_table.end(i); ++it) {
            CHECK(*it == debye_table.lookup(i, j++));
        }
    }
}