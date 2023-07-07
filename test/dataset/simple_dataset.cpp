#include <catch2/catch_test_macros.hpp>

#include <dataset/SimpleDataset.h>

struct fixture {
    SimpleDataset dataset = {{1, 2, 3}, {4, 5, 6}};
};

TEST_CASE("SimpleDataset::SimpleDataset") {
    SECTION("default constructor") {
        SimpleDataset dataset;
        CHECK(dataset.size() == 0);
        CHECK(dataset.size_rows() == 0);
        CHECK(dataset.size_cols() == 0);
    }

    SECTION("Dataset&") {
        SimpleDataset d1 = {{1, 2, 3}, {4, 5, 6}};
        SimpleDataset d2 = d1;
        CHECK(d1 == d2);
    }

    SECTION("unsigned int") {
        SimpleDataset dataset(10);
        CHECK(dataset.size() == 10);
        CHECK(dataset.size_rows() == 10);
        CHECK(dataset.size_cols() == 3);
    }

    SECTION("std::vector<double>, std::vector<double>, std::vector<double>") {
        SimpleDataset dataset({1, 2, 3}, {4, 5, 6}, {7, 8, 9});
        CHECK(dataset.size() == 3);
        CHECK(dataset.size_rows() == 3);
        CHECK(dataset.size_cols() == 3);
        CHECK(dataset.x() == std::vector{1, 2, 3});
        CHECK(dataset.y() == std::vector{4, 5, 6});
        CHECK(dataset.yerr() == std::vector{7, 8, 9});
    }

    SECTION("std::vector<double>, std::vector<double>") {
        SimpleDataset dataset({1, 2, 3}, {4, 5, 6});
        CHECK(dataset.size() == 3);
        CHECK(dataset.size_rows() == 3);
        CHECK(dataset.size_cols() == 3);
        CHECK(dataset.x() == std::vector{1, 2, 3});
        CHECK(dataset.y() == std::vector{4, 5, 6});
        CHECK(dataset.yerr() == std::vector{0, 0, 0});
    }

    SECTION("std::vector<double>, std::vector<double>, std::string, std::string") {
        SimpleDataset dataset({1, 2, 3}, {4, 5, 6}, "test1", "test2");
        CHECK(dataset.size() == 3);
        CHECK(dataset.size_rows() == 3);
        CHECK(dataset.size_cols() == 3);
        CHECK(dataset.x() == std::vector{1, 2, 3});
        CHECK(dataset.y() == std::vector{4, 5, 6});
        CHECK(dataset.yerr() == std::vector{0, 0, 0});
        CHECK(dataset.get_col_names()[0] == "test1");
        CHECK(dataset.get_col_names()[1] == "test2");
    }

    SECTION("io::ExistingFile") {
        SimpleDataset dataset("test/files/2epe.dat");
        CHECK(dataset.size() == 104);
        CHECK(dataset.size_rows() == 104);
        CHECK(dataset.size_cols() == 3);
    }
}

TEST_CASE_METHOD(fixture, "SimpleDataset::yerr") {    
    CHECK(dataset.yerr() == dataset.col(2));
}

TEST_CASE("SimpleDataset::load") {}

TEST_CASE("SimpleDataset::reduce") {}

TEST_CASE("SimpleDataset::operator=") {}

TEST_CASE("SimpleDataset::operator==") {}

// also tests span_x and span_y
TEST_CASE("SimpleDataset::get_limits") {
    SECTION("simple") {
        SimpleDataset dataset({1, 2, 3}, {4, 5, 6});
        CHECK(dataset.get_xlimits() == Limit(1, 3));
        CHECK(dataset.get_ylimits() == Limit(4, 6));
        CHECK(dataset.span_x() == dataset.get_xlimits());
        CHECK(dataset.span_y() == dataset.get_ylimits());
        CHECK(dataset.span_y_positive() == Limit(4, 6));
    }

    SECTION("complex") {
        SimpleDataset dataset({-2, 3, 5}, {-1, 4, 6});
        CHECK(dataset.get_xlimits() == Limit(-2, 5));
        CHECK(dataset.get_ylimits() == Limit(-1, 6));
        CHECK(dataset.span_x() == dataset.get_xlimits());
        CHECK(dataset.span_y() == dataset.get_ylimits());
        CHECK(dataset.span_y_positive() == Limit(0, 6));
    }
}

TEST_CASE("SimpleDataset::push_back") {
    SECTION("2D") {
        SimpleDataset dataset({1, 2, 3}, {4, 5, 6});
        dataset.push_back(4, 7);
        CHECK(dataset.size() == 4);
        CHECK(dataset.size_rows() == 4);
        CHECK(dataset.size_cols() == 3);
        CHECK(dataset.x() == std::vector{1, 2, 3, 4});
        CHECK(dataset.y() == std::vector{4, 5, 6, 7});
        CHECK(dataset.yerr() == std::vector{0, 0, 0, 0});
    } 

    SECTION("3D") {
        SimpleDataset dataset({1, 2, 3}, {4, 5, 6}, {7, 8, 9});
        dataset.push_back(4, 7, 10);
        CHECK(dataset.size() == 4);
        CHECK(dataset.size_rows() == 4);
        CHECK(dataset.size_cols() == 3);
        CHECK(dataset.x() == std::vector{1, 2, 3, 4});
        CHECK(dataset.y() == std::vector{4, 5, 6, 7});
        CHECK(dataset.yerr() == std::vector{7, 8, 9, 10});
    }
}

TEST_CASE("SimpleDataset::normalize") {}
TEST_CASE("SimpleDataset::scale_errors") {}
TEST_CASE("SimpleDataset::scale_y") {}
TEST_CASE("SimpleDataset::simulate_noise") {}
TEST_CASE("SimpleDataset::simulate_errors") {}
TEST_CASE("SimpleDataset::get_point") {}
TEST_CASE("SimpleDataset::find_minimum") {}
TEST_CASE("SimpleDataset::rebin") {}
TEST_CASE("SimpleDataset::generate_random_data") {}
TEST_CASE("SimpleDataset::mean") {}
TEST_CASE("SimpleDataset::weighted_mean") {}
TEST_CASE("SimpleDataset::std") {}
TEST_CASE("SimpleDataset::weighted_mean_error") {}
TEST_CASE("SimpleDataset::remove_consecutive_duplicates") {}