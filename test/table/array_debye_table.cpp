#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <table/ArrayDebyeTable.h>
#include <utility/Utility.h>
#include <constants/Constants.h>

auto& debye_table = table::ArrayDebyeTable::get_default_table();

TEST_CASE("DebyeTable::correct_values") {
    SECTION("default_table") {
        auto& q = constants::axes::q_axis;
        auto& d = constants::axes::d_axis;

        // prepare the q values for the intensity calculations
        REQUIRE(debye_table.size_d() == d.size());
        REQUIRE(debye_table.size_q() == q.size());
        for (unsigned int i = 0; i < q.size(); ++i) {
            for (unsigned int j = 0; j < d.size(); ++j) {
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