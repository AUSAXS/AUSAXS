#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <table/DebyeTable.h>
#include <utility/Axis.h>
#include <utility/Utility.h>
#include <settings/All.h>

TEST_CASE("DebyeTable::DebyeTable") {
    settings::general::verbose = false;
    SECTION("vector<double>&, vector<double>&") {
        std::vector<double> d = {1, 2, 3, 4, 5};
        std::vector<double> q = {1, 2, 3, 4, 5, 6};
        table::DebyeTable table(q, d);
        CHECK(table.size_d() == d.size());
        CHECK(table.size_q() == q.size());
    }
}

TEST_CASE("DebyeTable::lookup") {
    settings::general::verbose = false;
    std::vector<double> d = {20, 10, 30, 15, 10};
    std::vector<double> q = {0.01, 0.05, 0.1, 0.25, 0.5};
    table::DebyeTable table(q, d);

    auto func = [&] (double q, double d) {return std::sin(q*d)/(q*d);};

    CHECK_THAT(table.lookup(1, 1), Catch::Matchers::WithinAbs(func(q[1], d[1]), 1e-3));
    CHECK_THAT(table.lookup(2, 2), Catch::Matchers::WithinAbs(func(q[2], d[2]), 1e-3));
    CHECK_THAT(table.lookup(3, 3), Catch::Matchers::WithinAbs(func(q[3], d[3]), 1e-3));
    CHECK_THAT(table.lookup(4, 4), Catch::Matchers::WithinAbs(func(q[4], d[4]), 1e-3));
    CHECK_THAT(table.lookup(2, 0), Catch::Matchers::WithinAbs(func(q[2], d[0]), 1e-3));
    CHECK_THAT(table.lookup(4, 3), Catch::Matchers::WithinAbs(func(q[4], d[3]), 1e-3));
}

// size_d & size_q
TEST_CASE("DebyeTable::size") {
    settings::general::verbose = false;
    SECTION("(5, 5)") {
        std::vector<double> d = {1, 2, 3, 4, 5};
        std::vector<double> q = {1, 2, 3, 4, 5};
        table::DebyeTable table(q, d);
        CHECK(table.size_d() == d.size());
        CHECK(table.size_q() == q.size());
    }

    SECTION("(5, 7)") {
        std::vector<double> d = {1, 2, 3, 4, 5};
        std::vector<double> q = {1, 2, 3, 4, 5, 6, 7};
        table::DebyeTable table(q, d);
        CHECK(table.size_d() == d.size());
        CHECK(table.size_q() == q.size());
    }

    SECTION("empty") {
        std::vector<double> d = {};
        std::vector<double> q = {};
        table::DebyeTable table(q, d);
        CHECK(table.size_d() == d.size());
        CHECK(table.size_q() == q.size());
    }
}

TEST_CASE("DebyeTable::correct_values") {
    settings::general::verbose = false;
    settings::axes::qmax = 1.001;
    SECTION("default_table") {
        table::DebyeTable::reset_default_table();
        double width = settings::axes::distance_bin_width;
        auto q = Axis(settings::axes::qmin, settings::axes::qmax, settings::axes::bins).as_vector();
        std::vector<double> d = Axis(0, settings::axes::max_distance, settings::axes::max_distance/settings::axes::distance_bin_width).as_vector(0.5);
        d[0] = 0;

        // prepare the q values for the intensity calculations
        auto& table = table::DebyeTable::get_default_table();
        REQUIRE(table.size_d() == d.size());
        REQUIRE(table.size_q() == q.size());
        for (unsigned int i = 0; i < q.size(); ++i) {
            for (unsigned int j = 0; j < d.size(); ++j) {
                double qd = q[i]*d[j];
                if (qd < 1e-3) {
                    CHECK_THAT(table.lookup(i, j), Catch::Matchers::WithinAbs(1 - qd*qd/6 + qd*qd*qd*qd/120, 1e-3)); 
                    continue;
                }
                CHECK_THAT(table.lookup(i, j), Catch::Matchers::WithinAbs(std::sin(qd)/qd, 1e-3));
            }
        }
    }

    SECTION("non_default") {
        std::vector<double> q = {0.10, 0.12, 0.14, 0.16, 0.18, 0.20};
        std::vector<double> d = {0, 2, 4, 6, 8, 10};
        table::DebyeTable table(q, d);

        for (unsigned int i = 0; i < q.size(); ++i) {
            for (unsigned int j = 0; j < d.size(); ++j) {
                double qd = q[i]*d[j];
                if (qd < 1e-3) {
                    CHECK(table.lookup(i, j) == 1 - qd*qd/6 + qd*qd*qd*qd/120); 
                    continue;
                }
                CHECK_THAT(table.lookup(i, j), Catch::Matchers::WithinAbs(std::sin(qd)/qd, 1e-3));
            }
        }
    }
    settings::axes::qmax = 0.5;
}