#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

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

TEST_CASE("SimpleDataset::load") {
    SimpleDataset dataset(
        std::vector<double>{0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10}, 
        std::vector<double>{1,    2,    3,    4,    5,    6,    7,    8,    9,    10},
        std::vector<double>{1,    1,    1,    1,    1,    1,    1,    1,    1,    1}
    );
    std::string path = "temp/dataset/save.dat";
    dataset.save(path);

    SimpleDataset loaded_dataset(path);
    CHECK(dataset == loaded_dataset);
}

TEST_CASE("SimpleDataset::reduce") {
    SECTION("linear") {
        SimpleDataset dataset(
            std::vector<double>{0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11}, 
            std::vector<double>{1,    2,    3,    4,    5,    6,    7,    8,    9,    10,   11},
            std::vector<double>{1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1}
        );

        dataset.reduce(6);
        CHECK(dataset.size_rows() == 6);
        CHECK(dataset.size_cols() == 3);
        CHECK(dataset.x() == std::vector{0.01, 0.03, 0.05, 0.07, 0.09, 0.11});
    }

    SECTION("log") {
        SimpleDataset dataset(
            std::vector<double>{1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5}, 
            std::vector<double>{1,    2,    3,    4,    5,   6,   7,   8,   9,   10},
            std::vector<double>{1,    1,    1,    1,    1,   1,   1,   1,   1,   1}
        );

        dataset.reduce(6, true);
        CHECK(dataset.size_rows() == 6);
        CHECK(dataset.size_cols() == 3);
        CHECK(dataset.x() == std::vector{1, 2, 3, 4, 5, 6});
    }
}

TEST_CASE("SimpleDataset::operator=") {
    SimpleDataset dataset1({1, 2, 3}, {4, 5, 6});
    SimpleDataset dataset2({7, 8, 9}, {10, 11, 12});
    dataset1 = dataset2;
    CHECK(dataset1 == dataset2);
}

TEST_CASE("SimpleDataset::operator==") {
    SimpleDataset dataset1({1, 2, 3}, {4, 5, 6});
    SimpleDataset dataset2({7, 8, 9}, {10, 11, 12});
    CHECK(dataset1 != dataset2);

    dataset1 = dataset2;
    CHECK(dataset1 == dataset2);
}

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

    SECTION("empty") {
        SimpleDataset data;
        data.push_back(1, 2);
        data.push_back(3, 4);
        data.push_back(5, 6);
        data.push_back(Point2D(7, 8));

        CHECK(data.x() == Vector{1, 3, 5, 7});
        CHECK(data.y() == Vector{2, 4, 6, 8});
    }
}

TEST_CASE("SimpleDataset::normalize") {
    SimpleDataset dataset1(
        std::vector<double>{1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5}, 
        std::vector<double>{1,    2,    3,    4,    5,   6,   7,   8,   9,   10},
        std::vector<double>{1,    1,    1,    1,    1,   1,   1,   1,   1,   1}
    );
    auto dataset2 = dataset1;

    dataset2.normalize(1);
    CHECK(dataset2.y() == dataset1.y()*1e4);
}

TEST_CASE("SimpleDataset::scale_errors") {
    SimpleDataset dataset1(
        std::vector<double>{1e-4, 1e-3, 1e-2, 1e-1, 1e0}, 
        std::vector<double>{1,    2,    3,    4,    5},
        std::vector<double>{1,    1,    1,    1,    1}
    );
    auto dataset2 = dataset1;

    dataset2.scale_errors(2);
    CHECK(dataset2.yerr() == dataset1.yerr()*2);
}

TEST_CASE("SimpleDataset::scale_y") {
    SimpleDataset dataset1(
        std::vector<double>{1e-4, 1e-3, 1e-2, 1e-1, 1e0}, 
        std::vector<double>{1,    2,    3,    4,    5},
        std::vector<double>{1,    1,    1,    1,    1}
    );
    auto dataset2 = dataset1;

    dataset2.scale_y(2);
    CHECK(dataset2.y() == dataset1.y()*2);
}

TEST_CASE("SimpleDataset::simulate_noise") {
    SimpleDataset dataset1(
        std::vector<double>{1, 2, 3, 4, 5}, 
        std::vector<double>{1, 2, 3, 4, 5},
        std::vector<double>{1, 1, 1, 1, 1}
    );
    auto dataset2 = dataset1;

    std::vector<std::vector<double>> diff;
    for (unsigned int i = 0; i < 100; ++i) {
        dataset2.simulate_noise();
        diff.push_back(dataset2.y() - dataset1.y());
        CHECK_THAT(dataset1.mean(), Catch::Matchers::WithinAbs(dataset2.mean(), 1e-3));
        CHECK_THAT(dataset1.std(), Catch::Matchers::WithinAbs(dataset2.std(), 1e-3));
    }

    // Flatten the diff vector
    std::vector<double> flattened_diff;
    for (const auto& d : diff) {
        flattened_diff.insert(flattened_diff.end(), d.begin(), d.end());
    }

    // Calculate mean and standard deviation of the flattened_diff vector
    double mean = std::accumulate(flattened_diff.begin(), flattened_diff.end(), 0.0) / flattened_diff.size();
    double variance = 0.0;
    for (const auto& val : flattened_diff) {
        double diff = val - mean;
        variance += diff * diff;
    }
    variance /= flattened_diff.size();
    double std_dev = std::sqrt(variance);

    // Sort the flattened_diff vector
    std::sort(flattened_diff.begin(), flattened_diff.end());

    // Calculate the histogram
    const int num_bins = 20;
    std::vector<int> histogram(num_bins, 0);
    double min_value = flattened_diff.front();
    double max_value = flattened_diff.back();
    double bin_width = (max_value - min_value) / num_bins;

    for (const auto& val : flattened_diff) {
        int bin_index = static_cast<int>((val - min_value) / bin_width);
        histogram[bin_index]++;
    }

    // Check if the histogram is roughly symmetric
    bool is_symmetric = true;
    for (size_t i = 0; i < histogram.size() / 2; ++i) {
        if (histogram[i] != histogram[histogram.size() - i - 1]) {
            is_symmetric = false;
            break;
        }
    }

    // Check if the histogram follows a bell-shaped curve
    bool is_bell_shaped = true;
    for (size_t i = 1; i < histogram.size() - 1; ++i) {
        if (histogram[i] < histogram[i - 1] || histogram[i] < histogram[i + 1]) {
            is_bell_shaped = false;
            break;
        }
    }

    // Perform the checks
    CHECK(std_dev < 0.5);  // Check standard deviation is small
    CHECK(is_symmetric);  // Check histogram is roughly symmetric
    CHECK(is_bell_shaped);  // Check histogram follows a bell-shaped curve
}

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