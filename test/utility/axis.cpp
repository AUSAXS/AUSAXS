#include <catch2/catch_test_macros.hpp>

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
        CHECK(axis.get_bin(10) == 2);
        CHECK(axis.get_bin(11) == 2);
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
        CHECK(axis.get_bin(6) == 19);
    }
}

TEST_CASE("Axis::sub_axis") {
    SECTION("simple") {
        Axis axis(1, 10, 10);
        Axis sub_axis = axis.sub_axis(2, 8);
        CHECK(sub_axis.min == 2);
        CHECK(sub_axis.max == 8);
        CHECK(sub_axis.bins == 6);

        sub_axis = axis.sub_axis(0, 11);
        CHECK(sub_axis.min == 1);
        CHECK(sub_axis.max == 10);
        CHECK(sub_axis.bins == 10);

        sub_axis = axis.sub_axis(0.5, 5.5);
        CHECK(sub_axis.min == 1);
        CHECK(sub_axis.max == 5);
        CHECK(sub_axis.bins == 4);
    }

    SECTION("q_conversion") {
        auto q_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax).as_vector();
        CHECK(q_axis.size() == double(constants::axes::q_axis.bins)/2);
        CHECK(q_axis.front() == settings::axes::qmin);
        CHECK(q_axis.back() == settings::axes::qmax);

        for (unsigned int i = 0; i < q_axis.size(); ++i) {
            CHECK(q_axis[i] == constants::axes::q_axis.min + i*constants::axes::q_axis.width());
        }
    }
}