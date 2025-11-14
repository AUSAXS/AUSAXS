#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <table/DebyeTableManager.h>
#include <table/ArrayDebyeTable.h>
#include <table/VectorDebyeTable.h>
#include <constants/Constants.h>
#include <cmath>

using namespace ausaxs;
using namespace ausaxs::table;

static double sinc_approx(double q, double d) {
    double qd = q * d;
    if (std::abs(qd) < 1e-6) {
        return 1.0 - qd*qd/6.0 + qd*qd*qd*qd/120.0;
    }
    return std::sin(qd) / qd;
}

TEST_CASE("DebyeTableManager: default table is the shared default") {
    DebyeTableManager mgr;
    auto table_ptr = mgr.get_sinc_table();
    const auto &default_table = ArrayDebyeTable::get_default_table();

    REQUIRE(table_ptr != nullptr);
    CHECK(table_ptr == &default_table);
}

TEST_CASE("DebyeTableManager: custom axes produce expected sinc values") {
    DebyeTableManager mgr;

    const std::vector<double> q_axis{0.0, 0.1, 0.5};
    const std::vector<double> d_axis{0.0, 1.0, 2.0, 3.0};

    SECTION("accepts rvalue inputs (moves)") {
        mgr.set_q_axis(std::vector<double>(q_axis));
        mgr.set_d_axis(std::vector<double>(d_axis));

        auto tbl = mgr.get_sinc_table();
        REQUIRE(tbl != nullptr);
        CHECK(tbl != &ArrayDebyeTable::get_default_table());
        CHECK(tbl->size_q() == q_axis.size());
        CHECK(tbl->size_d() == d_axis.size());

        for (std::size_t i = 0; i < q_axis.size(); ++i) {
            for (std::size_t j = 0; j < d_axis.size(); ++j) {
                double expected = sinc_approx(q_axis[i], d_axis[j]);
                CHECK_THAT(tbl->lookup(static_cast<int>(i), static_cast<int>(j)), Catch::Matchers::WithinAbs(expected, 1e-6));
            }
        }
    }

    SECTION("accepts const-lvalue references (copies)") {
        mgr.set_q_axis(q_axis);
        mgr.set_d_axis(d_axis);

        auto tbl = mgr.get_sinc_table();
        REQUIRE(tbl != nullptr);
        CHECK(tbl != &ArrayDebyeTable::get_default_table());
        CHECK(tbl->size_q() == q_axis.size());
        CHECK(tbl->size_d() == d_axis.size());

        for (std::size_t i = 0; i < q_axis.size(); ++i) {
            for (std::size_t j = 0; j < d_axis.size(); ++j) {
                double expected = sinc_approx(q_axis[i], d_axis[j]);
                CHECK_THAT(tbl->lookup(static_cast<int>(i), static_cast<int>(j)), Catch::Matchers::WithinAbs(expected, 1e-6));
            }
        }
    }
}
