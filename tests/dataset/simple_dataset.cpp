#include "settings/GeneralSettings.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <dataset/SimpleDataset.h>
#include <math/Statistics.h>

struct fixture {
    SimpleDataset dataset = {{1, 2, 3}, {4, 5, 6}};
};

TEST_CASE("SimpleDataset::SimpleDataset") {
    settings::general::verbose = false;
    SECTION("default constructor") {
        SimpleDataset dataset;
        CHECK(dataset.size() == 0);
        CHECK(dataset.size_rows() == 0);
        CHECK(dataset.size_cols() == 3);
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

    SECTION("vector<double>, vector<double>, string, string") {
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
    settings::general::verbose = false;
    SECTION("simple") {
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

    SECTION("complex") {
        SimpleDataset data("test/files/2epe.dat");
        auto x = data.x();
        auto y = data.y();
        auto yerr = data.yerr();

        std::vector<double> validate_x = {9.81300045E-03, 1.06309997E-02, 1.14489999E-02, 1.22659998E-02, 1.30840000E-02, 1.39020002E-02, 1.47200003E-02, 1.55379996E-02, 1.63550004E-02, 1.71729997E-02};
        std::vector<double> validate_y = {6.67934353E-03, 7.27293547E-03, 8.74083303E-03, 9.22449585E-03, 9.13867634E-03, 9.21153929E-03, 9.37998667E-03, 8.67372658E-03, 9.23649967E-03, 9.22480784E-03};
        std::vector<double> validate_yerr = {1.33646582E-03, 1.01892441E-03, 8.62116576E-04, 7.71059655E-04, 6.87870081E-04, 6.30189374E-04, 4.98525158E-04, 4.69041377E-04, 4.46073769E-04, 4.26004088E-04};

        REQUIRE(x.size() == 104);
        REQUIRE(y.size() == 104);
        REQUIRE(yerr.size() == 104);
        for (unsigned int i = 0; i < validate_x.size(); i++) {
            CHECK_THAT(x[i], Catch::Matchers::WithinRel(validate_x[i]));
            CHECK_THAT(y[i], Catch::Matchers::WithinRel(validate_y[i]));
            CHECK_THAT(yerr[i], Catch::Matchers::WithinRel(validate_yerr[i]));
        }
    }
}

TEST_CASE("SimpleDataset::save") {
    settings::general::verbose = false;
    SECTION("same contents") {
        SimpleDataset data("test/files/2epe.dat");
        auto data2 = data;
        data.save("temp/dataset/2epe.dat");
        REQUIRE(data.size() == data2.size());
        REQUIRE(data.x().size() == data2.x().size());
        REQUIRE(data.y().size() == data2.y().size());
        REQUIRE(data.yerr().size() == data2.yerr().size());

        for (unsigned int i = 0; i < data.size(); i++) {
            CHECK_THAT(data.x(i), Catch::Matchers::WithinRel(data2.x(i), 1e-3));
            CHECK_THAT(data.y(i), Catch::Matchers::WithinRel(data2.y(i), 1e-3));
            CHECK_THAT(data.yerr(i), Catch::Matchers::WithinRel(data2.yerr(i), 1e-3));
        }
    }

    SECTION("accuracy") {
        SimpleDataset data;
        for (double x = 1.347e-01; x < 1.351e-01; x += 1e-6) {
            data.push_back(x, sin(x));
        }
        data.save("temp/dataset_io_accuracy.dat");

        SimpleDataset data2("temp/dataset_io_accuracy.dat");
        REQUIRE(data.size() == data2.size());
        for (unsigned int i = 0; i < data.size(); i++) {
            CHECK_THAT(data.x(i), Catch::Matchers::WithinAbs(data2.x(i), 1e-6));
            CHECK_THAT(data.y(i), Catch::Matchers::WithinAbs(data2.y(i), 1e-6));
        }
    }
}

TEST_CASE("SimpleDataset::reduce") {
    // we just want to check that it works as expected with perfectly linear data
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

    // we just want to check that it works as expected with perfectly logarithmic data
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
        CHECK(dataset.span_y_positive() == Limit(4, 6));
    }

    SECTION("old range tests") {
        std::vector<double> x = {1, 2, 3, 4, 5};
        std::vector<double> y = {10, 20, 30, 40, 50};
        SimpleDataset data(x, y);

        SECTION("limit") {
            SECTION("x") {
                data.limit_x(Limit(2, 3));
                CHECK(data.x() == std::vector<double>{2, 3});
                CHECK(data.y() == std::vector<double>{20, 30});
            }

            SECTION("y") {
                data.limit_y(Limit(20, 40));
                CHECK(data.x() == std::vector<double>{2, 3, 4});
                CHECK(data.y() == std::vector<double>{20, 30, 40});
            }
        }

        SECTION("spans") {
            SECTION("x") {
                auto span = data.span_x();
                CHECK(span == Limit(1, 5));
                CHECK(span == data.get_xlimits());
            }

            SECTION("y") {
                auto span = data.span_y();
                CHECK(span == Limit(10, 50));
                CHECK(span == data.get_ylimits());
            }

            SECTION("positive y") {
                y = {-6, -2, 1, 5, 8};
                data = SimpleDataset(x, y);
                auto span = data.span_y_positive();
                CHECK(span == Limit(1, 8));

                data = SimpleDataset();
                span = data.span_y_positive();
                CHECK(span == Limit(0, 0));
            }
        }
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
        std::vector<double>{1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5},
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

// TEST_CASE("SimpleDataset::simulate_noise", "[manual]") {
//     SimpleDataset dataset1(
//         std::vector<double>{1, 2, 3, 4, 5}, 
//         std::vector<double>{1, 2, 3, 4, 5},
//         std::vector<double>{1, 1, 1, 1, 1}
//     );
//     auto dataset2 = dataset1;

//     std::vector<std::vector<double>> diff;
//     for (unsigned int i = 0; i < 100; ++i) {
//         dataset2.simulate_noise();
//         diff.push_back(dataset2.y() - dataset1.y());
//     }

//     // Flatten the diff vector
//     std::vector<double> flattened_diff;
//     for (const auto& d : diff) {
//         flattened_diff.insert(flattened_diff.end(), d.begin(), d.end());
//     }

//     // Calculate mean and standard deviation of the flattened_diff vector
//     double mean = std::accumulate(flattened_diff.begin(), flattened_diff.end(), 0.0) / flattened_diff.size();
//     double variance = 0.0;
//     for (const auto& val : flattened_diff) {
//         double diff = val - mean;
//         variance += diff * diff;
//     }
//     variance /= flattened_diff.size();
//     double std_dev = std::sqrt(variance);

//     // Sort the flattened_diff vector
//     std::sort(flattened_diff.begin(), flattened_diff.end());

//     // Calculate the histogram
//     const int num_bins = 20;
//     std::vector<int> histogram(num_bins, 0);
//     double min_value = flattened_diff.front();
//     double max_value = flattened_diff.back();
//     double bin_width = (max_value - min_value) / num_bins;

//     for (const auto& val : flattened_diff) {
//         int bin_index = static_cast<int>((val - min_value) / bin_width);
//         histogram[bin_index]++;
//     }

//     // Check if the histogram is roughly symmetric
//     bool is_symmetric = true;
//     for (unsigned i = 0; i < histogram.size() / 2; ++i) {
//         if (histogram[i] != histogram[histogram.size() - i - 1]) {
//             is_symmetric = false;
//             break;
//         }
//     }

//     // Check if the histogram follows a bell-shaped curve
//     bool is_bell_shaped = true;
//     for (size_t i = 1; i < histogram.size() - 1; ++i) {
//         if (histogram[i] < histogram[i - 1] || histogram[i] < histogram[i + 1]) {
//             is_bell_shaped = false;
//             break;
//         }
//     }

//     // Perform the checks
//     CHECK(std_dev < 0.5);  // Check standard deviation is small
//     CHECK(is_symmetric);  // Check histogram is roughly symmetric
//     CHECK(is_bell_shaped);  // Check histogram follows a bell-shaped curve
// }

TEST_CASE("SimpleDataset::simulate_errors") {
    SimpleDataset dataset1(
        std::vector<double>{1, 2, 3, 4, 5}, 
        std::vector<double>{1, 2, 3, 4, 5},
        std::vector<double>{1, 1, 1, 1, 1}
    );
    auto dataset2 = dataset1;
    dataset1.simulate_errors();
    CHECK(dataset1.x() == dataset2.x());
    CHECK(dataset1.y() == dataset2.y());
    CHECK(dataset1.yerr() != dataset2.yerr());
}

TEST_CASE_METHOD(fixture, "SimpleDataset::get_point") {
    CHECK(dataset.get_point(0) == Point2D(1, 4, 0));
    CHECK(dataset.get_point(1) == Point2D(2, 5, 0));
    CHECK(dataset.get_point(2) == Point2D(3, 6, 0));
}

TEST_CASE("SimpleDataset::find_minimum") {
    settings::general::verbose = false;
    SECTION("simple") {
        SimpleDataset dataset(
            std::vector<double>{1, 2, 3, 4, 5}, 
            std::vector<double>{1, 2, 3, 4, 5},
            std::vector<double>{1, 1, 1, 1, 1}
        );
        CHECK(dataset.find_minimum() == Point2D(1, 1, 1));
    }

    SECTION("real data") {
        SimpleDataset dataset("test/files/SASDJQ4.dat");
        CHECK(dataset.find_minimum() == Point2D(0.39893, -0.00298, 0.06117));
    }
}

TEST_CASE("SimpleDataset::rebin") {
    SimpleDataset dataset("test/files/SASDJQ4.dat");
    REQUIRE(dataset.size() == 844);
    dataset.rebin();
    CHECK(dataset.size() < 400);
}

TEST_CASE("SimpleDataset::generate_random_data") {
    double val = GENERATE(10, 25, 50);
    unsigned int size = GENERATE(10, 20);

    auto dataset = SimpleDataset::generate_random_data(size, val);
    REQUIRE(dataset.size() == size);
    for (unsigned int i = 0; i < dataset.size(); ++i) {
        CHECK(dataset.x(i) == i);
        CHECK(dataset.y(i) >= -val);
        CHECK(dataset.y(i) <= val);
        CHECK(dataset.yerr(i) >= -val*0.1);
        CHECK(dataset.yerr(i) <= val*0.1);
    }
}

TEST_CASE("SimpleDataset::mean") {
    auto size = GENERATE(10, 25, 50);
    auto val = GENERATE(10, 5, 2);
    auto dataset = SimpleDataset::generate_random_data(size, val);
    CHECK(dataset.mean() == stats::mean(dataset.y()));
}

TEST_CASE("SimpleDataset::weighted_mean") {
    auto size = GENERATE(10, 25, 50);
    auto val = GENERATE(10, 5, 2);
    auto dataset = SimpleDataset::generate_random_data(size, val);
    CHECK(dataset.weighted_mean() == stats::weighted_mean(dataset.y(), dataset.yerr()));
}

TEST_CASE("SimpleDataset::std") {
    auto size = GENERATE(10, 25, 50);
    auto val = GENERATE(10, 5, 2);
    auto dataset = SimpleDataset::generate_random_data(size, val);
    CHECK(dataset.std() == stats::std(dataset.y()));
}

TEST_CASE("SimpleDataset::weighted_mean_error") {
    auto size = GENERATE(10, 25, 50);
    auto val = GENERATE(10, 5, 2);
    auto dataset = SimpleDataset::generate_random_data(size, val);
    CHECK(dataset.weighted_mean_error() == stats::weighted_mean_error(dataset.yerr()));
}

TEST_CASE("SimpleDataset::remove_consecutive_duplicates") {
    SECTION("many duplicates") {
        SimpleDataset dataset(
            std::vector<double>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, 
            std::vector<double>{1, 1, 2, 2, 3, 3, 4, 4, 5, 5},
            std::vector<double>{1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
        );
        dataset.remove_consecutive_duplicates();
        CHECK(dataset.x() == std::vector<double>{1, 3, 5, 7, 9});
        CHECK(dataset.y() == std::vector<double>{1, 2, 3, 4, 5});
        CHECK(dataset.yerr() == std::vector<double>{1, 1, 1, 1, 1});
    }

    SECTION("few duplicates") {
        std::vector<double> x = {1, 2, 3, 4, 5, 5, 5, 6, 7, 8, 9};
        std::vector<double> y = {10, 20, 30, 40, 50, 50, 50, 60, 70, 80, 90};
        SimpleDataset data(x, y);

        data.remove_consecutive_duplicates();
        CHECK(data.size() == 9);
        CHECK(data.x() == std::vector<double>({1, 2, 3, 4, 5, 6, 7, 8, 9}));
        CHECK(data.y() == std::vector<double>({10, 20, 30, 40, 50, 60, 70, 80, 90}));
    }

    SECTION("no duplicates") {
        std::vector<double> x = {1, 2, 3, 4, 5, 6, 7, 8, 9};
        std::vector<double> y = {10, 20, 30, 40, 50, 60, 70, 80, 90};
        SimpleDataset data(x, y);

        data.remove_consecutive_duplicates();
        CHECK(data.size() == 9);
        CHECK(data.x() == x);
        CHECK(data.y() == y);
    }
}