#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <dataset/SimpleDataset.h>
#include <dataset/PointSet.h>
#include <utility/Limit.h>

using namespace ausaxs;

TEST_CASE("SimpleDataset::SimpleDataset") {
    SECTION("default constructor") {
        SimpleDataset dataset;
        CHECK(dataset.size() == 0);
        CHECK(dataset.size_rows() == 0);
        CHECK(dataset.size_cols() == 3);
    }

    SECTION("copy constructor") {
        SimpleDataset d1 = {{1, 2, 3}, {4, 5, 6}};
        SimpleDataset d2 = d1;
        CHECK(d1 == d2);
    }

    SECTION("move constructor") {
        SimpleDataset d1 = {{1, 2, 3}, {4, 5, 6}};
        SimpleDataset d2 = std::move(d1);
        CHECK(d2.size() == 3);
        CHECK(d2.x() == std::vector{1, 2, 3});
        CHECK(d2.y() == std::vector{4, 5, 6});
    }

    SECTION("unsigned int") {
        SimpleDataset dataset(10);
        CHECK(dataset.size() == 10);
        CHECK(dataset.size_rows() == 10);
        CHECK(dataset.size_cols() == 3);
    }

    SECTION("vector<double>, vector<double>, vector<double>") {
        SimpleDataset dataset({1, 2, 3}, {4, 5, 6}, {7, 8, 9});
        CHECK(dataset.size() == 3);
        CHECK(dataset.size_rows() == 3);
        CHECK(dataset.size_cols() == 3);
        CHECK(dataset.x() == std::vector{1, 2, 3});
        CHECK(dataset.y() == std::vector{4, 5, 6});
        CHECK(dataset.yerr() == std::vector{7, 8, 9});
    }

    SECTION("vector<double>, vector<double>") {
        SimpleDataset dataset({1, 2, 3}, {4, 5, 6});
        CHECK(dataset.size() == 3);
        CHECK(dataset.size_rows() == 3);
        CHECK(dataset.size_cols() == 3);
        CHECK(dataset.x() == std::vector{1, 2, 3});
        CHECK(dataset.y() == std::vector{4, 5, 6});
        CHECK(dataset.yerr() == std::vector{0, 0, 0});
    }
}

TEST_CASE("SimpleDataset::yerr") {
    SimpleDataset dataset = {{1, 2, 3}, {4, 5, 6}};

    SECTION("column accessor") {
        CHECK(dataset.yerr() == dataset.col(2));
    }

    SECTION("element accessor") {
        SimpleDataset dataset2({1, 2, 3}, {4, 5, 6}, {7, 8, 9});
        CHECK(dataset2.yerr(0) == 7);
        CHECK(dataset2.yerr(1) == 8);
        CHECK(dataset2.yerr(2) == 9);
    }
}

TEST_CASE("SimpleDataset::push_back") {
    SimpleDataset dataset;

    SECTION("x, y, yerr") {
        dataset.push_back(1, 2, 3);
        CHECK(dataset.size() == 1);
        CHECK(dataset.x(0) == 1);
        CHECK(dataset.y(0) == 2);
        CHECK(dataset.yerr(0) == 3);

        dataset.push_back(4, 5, 6);
        CHECK(dataset.size() == 2);
        CHECK(dataset.x(1) == 4);
        CHECK(dataset.y(1) == 5);
        CHECK(dataset.yerr(1) == 6);
    }

    SECTION("Point2D") {
        Point2D point1(1, 2, 3);
        dataset.push_back(point1);
        CHECK(dataset.size() == 1);
        CHECK(dataset.x(0) == 1);
        CHECK(dataset.y(0) == 2);
        CHECK(dataset.yerr(0) == 3);

        Point2D point2(4, 5, 6);
        dataset.push_back(point2);
        CHECK(dataset.size() == 2);
        CHECK(dataset.x(1) == 4);
        CHECK(dataset.y(1) == 5);
        CHECK(dataset.yerr(1) == 6);
    }
}

TEST_CASE("SimpleDataset::get_point") {
    SimpleDataset dataset({1, 2, 3}, {4, 5, 6}, {7, 8, 9});

    CHECK(dataset.get_point(0) == Point2D(1, 4, 7));
    CHECK(dataset.get_point(1) == Point2D(2, 5, 8));
    CHECK(dataset.get_point(2) == Point2D(3, 6, 9));
}

TEST_CASE("SimpleDataset::span_x") {
    SECTION("empty dataset") {
        SimpleDataset dataset;
        Limit span = dataset.span_x();
        CHECK(span.min == 0);
        CHECK(span.max == 0);
    }

    SECTION("non-empty dataset") {
        SimpleDataset dataset({1, 5, 3, 9, 2}, {0, 0, 0, 0, 0});
        Limit span = dataset.span_x();
        CHECK(span.min == 1);
        CHECK(span.max == 9);
    }
}

TEST_CASE("SimpleDataset::span_y") {
    SECTION("empty dataset") {
        SimpleDataset dataset;
        Limit span = dataset.span_y();
        CHECK(span.min == 0);
        CHECK(span.max == 0);
    }

    SECTION("non-empty dataset") {
        SimpleDataset dataset({0, 0, 0, 0, 0}, {1, 5, 3, 9, 2});
        Limit span = dataset.span_y();
        CHECK(span.min == 1);
        CHECK(span.max == 9);
    }
}

TEST_CASE("SimpleDataset::span_y_positive") {
    SECTION("empty dataset") {
        SimpleDataset dataset;
        Limit span = dataset.span_y_positive();
        CHECK(span.min == 0);
        CHECK(span.max == 0);
    }

    SECTION("all positive") {
        SimpleDataset dataset({0, 0, 0, 0}, {1, 2, 3, 4});
        Limit span = dataset.span_y_positive();
        CHECK(span.min == 1);
        CHECK(span.max == 4);
    }

    SECTION("mixed positive and negative") {
        SimpleDataset dataset({0, 0, 0, 0, 0}, {-1, 2, -3, 4, 5});
        Limit span = dataset.span_y_positive();
        CHECK(span.min == 2);
        CHECK(span.max == 5);
    }
}

TEST_CASE("SimpleDataset::normalize") {
    SimpleDataset dataset({1, 2, 3}, {10, 20, 30});
    
    double factor = dataset.normalize(5);
    CHECK_THAT(factor, Catch::Matchers::WithinAbs(0.5, 1e-6));
    CHECK_THAT(dataset.y(0), Catch::Matchers::WithinAbs(5, 1e-6));
    CHECK_THAT(dataset.y(1), Catch::Matchers::WithinAbs(10, 1e-6));
    CHECK_THAT(dataset.y(2), Catch::Matchers::WithinAbs(15, 1e-6));
}

TEST_CASE("SimpleDataset::scale_y") {
    SimpleDataset dataset({1, 2, 3}, {10, 20, 30}, {1, 2, 3});
    
    dataset.scale_y(2);
    CHECK(dataset.y(0) == 20);
    CHECK(dataset.y(1) == 40);
    CHECK(dataset.y(2) == 60);
    CHECK(dataset.yerr(0) == 2);
    CHECK(dataset.yerr(1) == 4);
    CHECK(dataset.yerr(2) == 6);
}

TEST_CASE("SimpleDataset::scale_errors") {
    SimpleDataset dataset({1, 2, 3}, {10, 20, 30}, {1, 2, 3});
    
    dataset.scale_errors(3);
    CHECK(dataset.y(0) == 10);
    CHECK(dataset.y(1) == 20);
    CHECK(dataset.y(2) == 30);
    CHECK(dataset.yerr(0) == 3);
    CHECK(dataset.yerr(1) == 6);
    CHECK(dataset.yerr(2) == 9);
}

TEST_CASE("SimpleDataset::mean") {
    SimpleDataset dataset({1, 2, 3}, {10, 20, 30});
    CHECK_THAT(dataset.mean(), Catch::Matchers::WithinAbs(20, 1e-6));
}

TEST_CASE("SimpleDataset::std") {
    SimpleDataset dataset({1, 2, 3}, {10, 20, 30});
    double expected_std = 10.0; // sample std with ddof=1
    CHECK_THAT(dataset.std(), Catch::Matchers::WithinAbs(expected_std, 1e-6));
}

TEST_CASE("SimpleDataset::weighted_mean") {
    SimpleDataset dataset({1, 2, 3}, {10, 20, 30}, {1, 2, 3});
    double expected = (10/1.0 + 20/4.0 + 30/9.0) / (1/1.0 + 1/4.0 + 1/9.0);
    CHECK_THAT(dataset.weighted_mean(), Catch::Matchers::WithinAbs(expected, 1e-6));
}

TEST_CASE("SimpleDataset::weighted_mean_error") {
    SimpleDataset dataset({1, 2, 3}, {10, 20, 30}, {1, 2, 3});
    double expected = std::sqrt(1.0 / (1/1.0 + 1/4.0 + 1/9.0));
    CHECK_THAT(dataset.weighted_mean_error(), Catch::Matchers::WithinAbs(expected, 1e-6));
}

TEST_CASE("SimpleDataset::operator=") {
    SimpleDataset dataset1({1, 2, 3}, {4, 5, 6});
    SimpleDataset dataset2({7, 8, 9}, {10, 11, 12});
    
    dataset1 = dataset2;
    CHECK(dataset1 == dataset2);
    CHECK(dataset1.x() == std::vector{7, 8, 9});
    CHECK(dataset1.y() == std::vector{10, 11, 12});
}

TEST_CASE("SimpleDataset::operator==") {
    SimpleDataset dataset1({1, 2, 3}, {4, 5, 6});
    SimpleDataset dataset2({7, 8, 9}, {10, 11, 12});
    
    CHECK(dataset1 != dataset2);
    
    dataset1 = dataset2;
    CHECK(dataset1 == dataset2);
}

TEST_CASE("SimpleDataset::find_minimum") {
    SECTION("simple") {
        SimpleDataset dataset({1, 2, 3, 4, 5}, {5, 2, 4, 1, 3}, {1, 1, 1, 1, 1});
        CHECK(dataset.find_minimum() == Point2D(4, 1, 1));
    }

    SECTION("first element") {
        SimpleDataset dataset({1, 2, 3}, {1, 2, 3}, {1, 1, 1});
        CHECK(dataset.find_minimum() == Point2D(1, 1, 1));
    }

    SECTION("last element") {
        SimpleDataset dataset({1, 2, 3}, {3, 2, 1}, {1, 1, 1});
        CHECK(dataset.find_minimum() == Point2D(3, 1, 1));
    }
}

TEST_CASE("SimpleDataset::generate_random_data") {
    SECTION("two arguments") {
        SimpleDataset dataset = SimpleDataset::generate_random_data(10, 5);
        CHECK(dataset.size() == 10);
        CHECK(dataset.size_cols() == 3);
        
        for (unsigned int i = 0; i < dataset.size(); i++) {
            CHECK(dataset.x(i) == i);
            CHECK(dataset.y(i) >= -5);
            CHECK(dataset.y(i) <= 5);
        }
    }

    SECTION("three arguments") {
        SimpleDataset dataset = SimpleDataset::generate_random_data(10, -3, 7);
        CHECK(dataset.size() == 10);
        CHECK(dataset.size_cols() == 3);
        
        for (unsigned int i = 0; i < dataset.size(); i++) {
            CHECK(dataset.x(i) == i);
            CHECK(dataset.y(i) >= -3);
            CHECK(dataset.y(i) <= 7);
        }
    }
}

TEST_CASE("SimpleDataset::remove_consecutive_duplicates") {
    SECTION("no duplicates") {
        SimpleDataset dataset({1, 2, 3}, {1, 2, 3});
        dataset.remove_consecutive_duplicates();
        CHECK(dataset.size() == 3);
    }

    SECTION("with duplicates") {
        SimpleDataset dataset({1, 2, 3, 4, 5, 6, 7}, {1, 1, 2, 2, 2, 3, 3});
        dataset.remove_consecutive_duplicates();
        CHECK(dataset.size() == 3);
        CHECK(dataset.x() == Vector{1, 3, 6});
        CHECK(dataset.y() == Vector{1, 2, 3});
    }

    SECTION("all same") {
        SimpleDataset dataset({1, 2, 3, 4}, {5, 5, 5, 5});
        dataset.remove_consecutive_duplicates();
        CHECK(dataset.size() == 1);
        CHECK(dataset.x(0) == 1);
        CHECK(dataset.y(0) == 5);
    }
}

TEST_CASE("SimpleDataset::reduce") {
    SECTION("linear reduction") {
        SimpleDataset dataset({1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, 
                            {10, 20, 30, 40, 50, 60, 70, 80, 90, 100});
        dataset.reduce(5, false);
        CHECK(dataset.size() <= 5);
    }

    SECTION("logarithmic reduction") {
        std::vector<double> x(100);
        std::vector<double> y(100);
        for (int i = 0; i < 100; i++) {
            x[i] = std::pow(10, i * 0.05);
            y[i] = i;
        }
        SimpleDataset dataset(x, y);
        dataset.reduce(20, true);
        CHECK(dataset.size() <= 20);
    }

    SECTION("target larger than size") {
        SimpleDataset dataset({1, 2, 3}, {4, 5, 6});
        dataset.reduce(10, false);
        CHECK(dataset.size() == 3);
    }

    SECTION("linear") {
        SimpleDataset dataset(
            std::vector<double>{0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11}, 
            std::vector<double>{1,    2,    3,    4,    5,    6,    7,    8,    9,    10,   11},
            std::vector<double>{1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1}
        );

        dataset.reduce(6);
        CHECK(dataset.size_rows() < 11);
        CHECK(dataset.size_cols() == 3);
        CHECK(dataset.x() == std::vector{0.01, 0.03, 0.05, 0.07, 0.09, 0.11});
    }

    SECTION("log") {
        SimpleDataset dataset(
            std::vector<double>{1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6}, 
            std::vector<double>{1,    2,    3,    4,    5,   6,   7,   8,   9,   10,  11},
            std::vector<double>{1,    1,    1,    1,    1,   1,   1,   1,   1,   1,   1}
        );

        dataset.reduce(6, true);
        CHECK(dataset.size_rows() == 6);
        CHECK(dataset.size_cols() == 3);
        CHECK(dataset.x() == std::vector{1e-4, 1e-2, 1e0, 1e2, 1e4, 1e6});
    }
}
