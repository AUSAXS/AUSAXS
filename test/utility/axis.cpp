#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <utility/Axis.h>
#include <utility/Limit.h>
#include <constants/Constants.h>
#include <settings/HistogramSettings.h>

TEST_CASE("Axis::Axis") {
    SECTION("default") {
        Axis axis;
        CHECK(axis.span() == 0);
        CHECK(axis.empty());
    }

    SECTION("Limit&, int") {
        Limit limit(1, 10);
        Axis axis(limit, 3);
        CHECK(axis.width() == 3);
        CHECK(axis.span() == 9);
        CHECK(axis.step() == 3);
        CHECK_FALSE(axis.empty());
    }

    SECTION("int, double, double") {
        Axis axis(1, 10, 3);
        CHECK(axis.width() == 3);
        CHECK(axis.span() == 9);
        CHECK(axis.step() == 3);
        CHECK_FALSE(axis.empty());
    }
}

TEST_CASE("Axis::operator=") {
    Axis axis1(1, 9, 3);
    Axis axis2;
    axis2 = axis1;
    CHECK(axis1 == axis2);
}

TEST_CASE("Axis::operator==") {
    Axis axis1(1, 9, 3);
    Axis axis2(1, 9, 3);
    CHECK(axis1 == axis2);

    Axis axis3(1, 9, 4);
    CHECK_FALSE(axis1 == axis3);
}

TEST_CASE("Axis::operator!=") {
    Axis axis1(1, 9, 3);
    Axis axis2(1, 9, 3);
    CHECK_FALSE(axis1 != axis2);

    Axis axis3(1, 9, 4);
    CHECK(axis1 != axis3);
}

TEST_CASE("Axis::width") {
    Axis axis(1, 10, 3);
    CHECK(axis.width() == 3);
}

TEST_CASE("Axis::span") {
    Axis axis(1, 9, 3);
    CHECK(axis.span() == 8);
}

TEST_CASE("Axis::step") {
    Axis axis(1, 10, 3);
    CHECK(axis.step() == 3);
}

TEST_CASE("Axis::resize") {
    Axis axis(1, 10, 3);
    axis.resize(2);
    CHECK(axis.width() == 3);
    CHECK(axis.span() == 6);
    CHECK(axis.max == 7);
}

TEST_CASE("Axis::as_vector") {
    SECTION("use_center_values = false") {
        Axis axis(1, 10, 3);
        std::vector<double> vec = axis.as_vector();
        CHECK(vec.size() == 3);
        CHECK(vec[0] == 1);
        CHECK(vec[1] == 4);
        CHECK(vec[2] == 7);
    }

    SECTION("use_center_values = true") {
        Axis axis(1, 10, 3);
        std::vector<double> vec = axis.as_vector(0.5);
        CHECK(vec.size() == 3);
        CHECK(vec[0] == 2.5);
        CHECK(vec[1] == 5.5);
        CHECK(vec[2] == 8.5);
    }
}

TEST_CASE("Axis::empty") {
    Axis axis;
    CHECK(axis.empty());

    axis = Axis(1, 10, 3);
    CHECK_FALSE(axis.empty());
}

TEST_CASE("Axis::limits") {
    Axis axis(1, 10, 3);
    CHECK(axis.limits().min == 1);
    CHECK(axis.limits().max == 10);
}

TEST_CASE("Axis::get_bin") {
    SECTION("empty") {
        Axis axis;
        CHECK(axis.get_bin(1) == 0);
    }

    SECTION("int") {
        Axis axis(1, 10, 3);
        CHECK(axis.get_bin(1) == 0);
        CHECK(axis.get_bin(2) == 0);
        CHECK(axis.get_bin(3) == 0);
        CHECK(axis.get_bin(4) == 1);
        CHECK(axis.get_bin(5) == 1);
        CHECK(axis.get_bin(6) == 1);
        CHECK(axis.get_bin(7) == 2);
        CHECK(axis.get_bin(8) == 2);
        CHECK(axis.get_bin(9) == 2);
        CHECK(axis.get_bin(10) == 3);
        CHECK(axis.get_bin(11) == 3);
    }

    SECTION("double") {
        Axis axis(1, 5, 20);
        for (unsigned int i = 0; i < axis.bins; ++i) {
            if (axis.get_bin(axis.min + i*axis.width()) != i) {
                std::cout << "i: " << i << std::endl;
                std::cout << "\tget_bin(" << axis.min + i*axis.width() << "): " << axis.get_bin(axis.min + i*axis.width()) << std::endl;
                CHECK(false);
            }
        }
        CHECK(axis.get_bin(0) == 0);
        CHECK(axis.get_bin(0.24) == 0);
        CHECK(axis.get_bin(6) == 20);
    }
}

TEST_CASE("Axis::get_bin_value") {
    SECTION("empty") {
        Axis axis;
        CHECK(axis.get_bin_value(0) == 0);
    }

    SECTION("simple") {
        Axis axis(0, 10, 10);
        for (unsigned int i = 0; i < axis.bins; ++i) {
            CHECK(axis.get_bin_value(i) == i);
        }        

        axis = Axis(1, 10, 3);
        CHECK(axis.get_bin_value(0) == 1);
        CHECK(axis.get_bin_value(1) == 4);
        CHECK(axis.get_bin_value(2) == 7);
        CHECK(axis.get_bin_value(3) == 10);
    }

    SECTION("complex") {
        Axis axis(1, 5, 20);
        for (unsigned int i = 0; i < axis.bins; ++i) {
            CHECK_THAT(axis.get_bin_value(i), Catch::Matchers::WithinAbs(axis.min + i*axis.width(), 1e-6));
        }
    }

    SECTION("d_axis") {
        auto& d_axis = constants::axes::d_axis;
        for (unsigned int i = 0; i < d_axis.bins; ++i) {
            CHECK_THAT(d_axis.get_bin_value(i), Catch::Matchers::WithinAbs(0.5*i, 1e-6));
        }
    }
}

TEST_CASE("Axis::sub_axis") {
    SECTION("simple") {
        Axis axis(1, 10, 9);
        Axis sub_axis = axis.sub_axis(2, 8);
        CHECK(sub_axis.min == 2);
        CHECK(sub_axis.max == 8);
        CHECK(sub_axis.bins == 6);

        sub_axis = axis.sub_axis(0, 11);
        CHECK(sub_axis.min == 1);
        CHECK(sub_axis.max == 10);
        CHECK(sub_axis.bins == 9);

        sub_axis = axis.sub_axis(0.5, 5.5);
        CHECK(sub_axis.min == 1);
        CHECK(sub_axis.max == 5);
        CHECK(sub_axis.bins == 4);
    }

    SECTION("q_conversion") {
        auto q_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);

        // settings::axes::qmax is not an exact entry in constants::axes::q_axis, so we can only check that we are close
        CHECK_THAT(q_axis.max, Catch::Matchers::WithinAbs(settings::axes::qmax, 1e-2)); 
        CHECK(q_axis.min == settings::axes::qmin);
        CHECK_THAT(q_axis.bins, Catch::Matchers::WithinAbs(constants::axes::q_axis.bins/2, 1.01)); 

        auto qvals = q_axis.as_vector();
        for (unsigned int i = 0; i < q_axis.bins; ++i) {
            CHECK(qvals[i] == constants::axes::q_vals[i]);
        }
    }
}