#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <dataset/Dataset.h>
#include <dataset/PointSet.h>
#include <dataset/SimpleDataset.h>
#include <settings/GeneralSettings.h>

#include <support/temp_file.h>

using namespace ausaxs;

TEST_CASE("SimpleDataset::SimpleDataset") {
    SECTION("default constructor") {
        SimpleDataset dataset;
        CHECK(dataset.empty());
        CHECK(dataset.size_rows() == 0);
        CHECK(dataset.size_cols() == 3);
    }

    SECTION("copy constructor") {
        SimpleDataset d1({1, 2, 3}, {4, 5, 6}, {1, 1, 1});
        SimpleDataset d2 = d1;
        CHECK(d1 == d2);
    }

    SECTION("move constructor") {
        SimpleDataset d1({1, 2, 3}, {4, 5, 6}, {1, 1, 1});
        SimpleDataset d2 = std::move(d1);
        CHECK(d2.size() == 3);
        CHECK(d2.x() == std::vector{1, 2, 3});
        CHECK(d2.y() == std::vector{4, 5, 6});
    }

    SECTION("int") {
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

    SECTION("Dataset") {
        SECTION("two columns") {
            Dataset d{{1, 2, 3}, {4, 5, 6}};
            CHECK_THROWS(SimpleDataset(d));
        }

        SECTION("four columns") {
            Dataset d{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}, {10, 11, 12}};
            SimpleDataset dataset(d);
            CHECK(dataset.size_cols() == 3);
            CHECK(dataset.x() == std::vector{1, 2, 3});
            CHECK(dataset.y() == std::vector{4, 5, 6});
            CHECK(dataset.yerr() == std::vector{7, 8, 9});
        }
    }

    SECTION("io::ExistingFile") {
        settings::general::verbose = false;
        SimpleDataset dataset("tests/files/2epe.dat");
        CHECK(dataset.size() == 104);
        CHECK(dataset.size_rows() == 104);
        CHECK(dataset.size_cols() == 3);
    }
}

TEST_CASE("SimpleDataset::yerr") {
    SimpleDataset dataset({1, 2, 3}, {4, 5, 6}, {7, 8, 9});

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

TEST_CASE("SimpleDataset::normalize") {
    SimpleDataset dataset({1, 2, 3}, {10, 20, 30}, {1, 1, 1});
    
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
    SimpleDataset dataset1({1, 2, 3}, {4, 5, 6}, {1, 1, 1});
    SimpleDataset dataset2({7, 8, 9}, {10, 11, 12}, {1, 1, 1});
    
    dataset1 = dataset2;
    CHECK(dataset1 == dataset2);
    CHECK(dataset1.x() == std::vector{7, 8, 9});
    CHECK(dataset1.y() == std::vector{10, 11, 12});
}

TEST_CASE("SimpleDataset::operator==") {
    SimpleDataset dataset1({1, 2, 3}, {4, 5, 6}, {1, 1, 1});
    SimpleDataset dataset2({7, 8, 9}, {10, 11, 12}, {1, 1, 1});
    
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
        
        for (int i = 0; i < dataset.size(); i++) {
            CHECK(dataset.x(i) == i);
            CHECK(dataset.y(i) >= -5);
            CHECK(dataset.y(i) <= 5);
        }
    }

    SECTION("three arguments") {
        SimpleDataset dataset = SimpleDataset::generate_random_data(10, -3, 7);
        CHECK(dataset.size() == 10);
        CHECK(dataset.size_cols() == 3);
        
        for (int i = 0; i < dataset.size(); i++) {
            CHECK(dataset.x(i) == i);
            CHECK(dataset.y(i) >= -3);
            CHECK(dataset.y(i) <= 7);
        }
    }
}

TEST_CASE("SimpleDataset::remove_consecutive_duplicates") {
    SECTION("no duplicates") {
        SimpleDataset dataset({1, 2, 3}, {1, 2, 3}, {1, 1, 1});
        dataset.remove_consecutive_duplicates();
        CHECK(dataset.size() == 3);
    }

    SECTION("with duplicates") {
        SimpleDataset dataset({1, 2, 3, 4, 5, 6, 7}, {1, 1, 2, 2, 2, 3, 3}, {1, 1, 1, 1, 1, 1, 1});
        dataset.remove_consecutive_duplicates();
        CHECK(dataset.size() == 3);
        CHECK(dataset.x() == Vector{1, 3, 6});
        CHECK(dataset.y() == Vector{1, 2, 3});
    }

    SECTION("all same") {
        SimpleDataset dataset({1, 2, 3, 4}, {5, 5, 5, 5}, {1, 1, 1, 1});
        dataset.remove_consecutive_duplicates();
        CHECK(dataset.size() == 1);
        CHECK(dataset.x(0) == 1);
        CHECK(dataset.y(0) == 5);
    }
}

TEST_CASE("SimpleDataset::reduce") {
    SECTION("linear reduction") {
        SimpleDataset dataset({1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, 
                            {10, 20, 30, 40, 50, 60, 70, 80, 90, 100},
                            {1, 1, 1, 1, 1, 1, 1, 1, 1, 1});
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
        SimpleDataset dataset(x, y, std::vector<double>(100, 1));
        dataset.reduce(20, true);
        CHECK(dataset.size() <= 20);
    }

    SECTION("target larger than size") {
        SimpleDataset dataset({1, 2, 3}, {4, 5, 6}, {1, 1, 1});
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

TEST_CASE("SimpleDataset::load") {
    settings::general::verbose = false;
    SECTION("simple") {
        SimpleDataset dataset(
            std::vector<double>{0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10}, 
            std::vector<double>{1,    2,    3,    4,    5,    6,    7,    8,    9,    10},
            std::vector<double>{1,    1,    1,    1,    1,    1,    1,    1,    1,    1}
        );
        test::TempFile path(".dat");
        dataset.save(path);

        SimpleDataset loaded_dataset(path);
        CHECK(dataset == loaded_dataset);
    }

    SECTION("complex") {
        SimpleDataset data("tests/files/2epe.dat");
        auto x = data.x();
        auto y = data.y();
        auto yerr = data.yerr();

        std::vector<double> validate_x = {9.81300045E-03, 1.06309997E-02, 1.14489999E-02, 1.22659998E-02, 1.30840000E-02, 1.39020002E-02, 1.47200003E-02, 1.55379996E-02, 1.63550004E-02, 1.71729997E-02};
        std::vector<double> validate_y = {6.67934353E-03, 7.27293547E-03, 8.74083303E-03, 9.22449585E-03, 9.13867634E-03, 9.21153929E-03, 9.37998667E-03, 8.67372658E-03, 9.23649967E-03, 9.22480784E-03};
        std::vector<double> validate_yerr = {1.33646582E-03, 1.01892441E-03, 8.62116576E-04, 7.71059655E-04, 6.87870081E-04, 6.30189374E-04, 4.98525158E-04, 4.69041377E-04, 4.46073769E-04, 4.26004088E-04};

        REQUIRE(x.size() == 104);
        REQUIRE(y.size() == 104);
        REQUIRE(yerr.size() == 104);
        for (int i = 0; i < static_cast<int>(validate_x.size()); i++) {
            CHECK_THAT(x[i], Catch::Matchers::WithinRel(validate_x[i]));
            CHECK_THAT(y[i], Catch::Matchers::WithinRel(validate_y[i]));
            CHECK_THAT(yerr[i], Catch::Matchers::WithinRel(validate_yerr[i]));
        }
    }
}

TEST_CASE("SimpleDataset::save") {
    settings::general::verbose = false;
    SECTION("same contents") {
        SimpleDataset data("tests/files/2epe.dat");
        auto data2 = data;
        data.save(test::TempFile(".dat"));
        REQUIRE(data.size() == data2.size());
        REQUIRE(data.x().size() == data2.x().size());
        REQUIRE(data.y().size() == data2.y().size());
        REQUIRE(data.yerr().size() == data2.yerr().size());

        for (int i = 0; i < data.size(); i++) {
            CHECK_THAT(data.x(i), Catch::Matchers::WithinRel(data2.x(i), 1e-3));
            CHECK_THAT(data.y(i), Catch::Matchers::WithinRel(data2.y(i), 1e-3));
            CHECK_THAT(data.yerr(i), Catch::Matchers::WithinRel(data2.yerr(i), 1e-3));
        }
    }

    SECTION("accuracy") {
        SimpleDataset data;
        for (double x = 1.347e-01; x < 1.351e-01; x += 1e-6) {
            data.push_back(x, sin(x), 1);
        }
        test::TempFile path(".dat");
        data.save(path);

        SimpleDataset data2(path);
        REQUIRE(data.size() == data2.size());
        for (int i = 0; i < data.size(); i++) {
            CHECK_THAT(data.x(i), Catch::Matchers::WithinAbs(data2.x(i), 1e-6));
            CHECK_THAT(data.y(i), Catch::Matchers::WithinAbs(data2.y(i), 1e-6));
        }
    }
}

TEST_CASE("SimpleDataset::simulate_errors") {
    SimpleDataset dataset1(
        std::vector<double>{1, 2, 3, 4, 5}, 
        std::vector<double>{1, 2, 3, 4, 5},
        std::vector<double>{0, 0, 0, 0, 0}
    );
    auto dataset2 = dataset1;

    dataset2.simulate_errors();
    CHECK_FALSE(dataset2.yerr() == dataset1.yerr());
}

TEST_CASE("SimpleDataset::rebin") {
    settings::general::verbose = false;
    SimpleDataset dataset("tests/files/SASDJQ4.dat");
    REQUIRE(dataset.size() == 844);
    dataset.rebin();
    CHECK(dataset.size() < 400);

    // below q = 0.03 the fold factor is 1, and folding a point into itself is the identity
    SECTION("fold factor 1 is the identity") {
        double err_scale = GENERATE(0., 0.01);
        std::vector<double> x, y, yerr;
        for (int i = 0; i < 20; ++i) {
            x.push_back(0.001*(i+1));
            y.push_back(i);
            yerr.push_back(err_scale*(i+1));
        }
        SimpleDataset unfolded(x, y, yerr);
        unfolded.rebin();

        REQUIRE(unfolded.size() == static_cast<int>(x.size()));
        for (int i = 0; i < static_cast<int>(x.size()); ++i) {
            CHECK_THAT(unfolded.x(i),    Catch::Matchers::WithinAbs(x[i],    1e-9));
            CHECK_THAT(unfolded.y(i),    Catch::Matchers::WithinAbs(y[i],    1e-9));
            CHECK_THAT(unfolded.yerr(i), Catch::Matchers::WithinAbs(yerr[i], 1e-9));
        }
    }

    // every input point must be folded into exactly one output point
    SECTION("no input points are dropped") {
        // a band with a constant fold factor, so the expected grouping is known exactly
        auto [qmin, fold] = GENERATE(
            std::make_pair(0.03, 2),
            std::make_pair(0.06, 4),
            std::make_pair(0.10, 8)
        );

        int N = 8*fold;
        std::vector<double> x, y, yerr(N, 1);
        for (int i = 0; i < N; ++i) {
            x.push_back(qmin + 1e-5*(i+1));
            y.push_back(i);                 // y_i = i, so each bin reports the mean of its indices
        }
        SimpleDataset folded(x, y, yerr);
        folded.rebin();

        REQUIRE(folded.size() == N/fold);
        for (int i = 0; i < folded.size(); ++i) {
            double first = i*fold;
            CHECK_THAT(folded.y(i), Catch::Matchers::WithinAbs(first + (fold-1)/2., 1e-9));
        }
    }
}
