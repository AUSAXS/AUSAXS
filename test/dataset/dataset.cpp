#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <dataset/Dataset.h>
#include <math/Matrix.h>
#include <dataset/PointSet.h>
#include <utility/Limit.h>
#include <io/ExistingFile.h>

#include <fstream>

struct fixture {
    std::vector<double> x = {   1,   2,   3,   4,   5,   6,   7,   8,   9};
    std::vector<double> y = {  -6,  -4,  -1,   2,   1,   3,   6,   7,   9};
    Dataset dataset = Dataset({x, y});
};

TEST_CASE("Dataset::Dataset") {
    SECTION("default") {
        Dataset dataset;
        CHECK(dataset.size() == 0);
        CHECK(dataset.size_rows() == 0);
        CHECK(dataset.size_cols() == 0);
    }

    SECTION("Dataset&") {
        std::vector<std::vector<double>> cols = {{1, 2, 3}, {4, 5, 6}};
        std::vector<std::string> col_names = {"a", "b"};
        Dataset dataset(cols, col_names);
        Dataset dataset2(dataset);
        CHECK(dataset == dataset2);
    }

    SECTION("Dataset&&") {
        std::vector<std::vector<double>> cols = {{1, 2, 3}, {4, 5, 6}};
        std::vector<std::string> col_names = {"a", "b"};
        Dataset dataset(cols, col_names);
        Dataset dataset2(std::move(dataset));
        REQUIRE(dataset2.size() == 3);
        REQUIRE(dataset2.size_rows() == 3);
        REQUIRE(dataset2.size_cols() == 2);
        CHECK(dataset2.get_col_names() == col_names);
        CHECK(dataset2.col("a") == Vector<double>{1, 2, 3});
        CHECK(dataset2.col("b") == Vector<double>{4, 5, 6});
    }

    SECTION("Matrix&&") {
        Matrix<double> m(2, 2);
        Dataset dataset(std::move(m));
        CHECK(dataset.size() == 2);
        CHECK(dataset.size_rows() == 2);
        CHECK(dataset.size_cols() == 2);
    }

    SECTION("vector<string>&") {
        std::vector<std::string> col_names = {"a", "b", "c"};
        Dataset dataset(col_names);
        CHECK(dataset.size() == 0);
        CHECK(dataset.size_rows() == 0);
        CHECK(dataset.size_cols() == 3);
        CHECK(dataset.get_col_names() == col_names);
    }

    SECTION("vector<vector<double>>&, vector<string>&") {
        std::vector<std::vector<double>> cols = {{1, 2, 3}, {4, 5, 6}};
        std::vector<std::string> col_names = {"a", "b"};
        Dataset dataset(cols, col_names);
        CHECK(dataset.size() == 3);
        CHECK(dataset.size_rows() == 3);
        CHECK(dataset.size_cols() == 2);
        CHECK(dataset.get_col_names() == col_names);
        CHECK(dataset.col("a") == Vector<double>{1, 2, 3});
        CHECK(dataset.col("b") == Vector<double>{4, 5, 6});
    }

    SECTION("vector<vector<double>>&") {
        std::vector<std::vector<double>> cols = {{1, 2, 3}, {4, 5, 6}};
        Dataset dataset(cols);
        CHECK(dataset.size() == 3);
        CHECK(dataset.size_rows() == 3);
        CHECK(dataset.size_cols() == 2);
        CHECK(dataset.col(0) == Vector<double>{1, 2, 3});
        CHECK(dataset.col(1) == Vector<double>{4, 5, 6});
    }

    SECTION("unsigned int, unsigned int") {
        Dataset dataset(2, 3);
        CHECK(dataset.size() == 2);
        CHECK(dataset.size_rows() == 2);
        CHECK(dataset.size_cols() == 3);
    }

    SECTION("ExistingFile&") {
        io::ExistingFile file("test/files/2epe.dat");
        Dataset dataset(file);
        CHECK(dataset.size() == 104);
        CHECK(dataset.size_rows() == 104);
        CHECK(dataset.size_cols() == 4);
    }
}

TEST_CASE_METHOD(fixture, "Dataset::col") {
    CHECK(dataset.col(0) == Vector<double>{1, 2, 3, 4, 5, 6, 7, 8, 9});
    CHECK(dataset.col(1) == Vector<double>{-6, -4, -1, 2, 1, 3, 6, 7, 9});
}

TEST_CASE_METHOD(fixture, "Dataset::row") {
    CHECK(dataset.row(0) == Vector<double>{1, -6});
    CHECK(dataset.row(1) == Vector<double>{2, -4});
    CHECK(dataset.row(2) == Vector<double>{3, -1});
    CHECK(dataset.row(3) == Vector<double>{4, 2});
    CHECK(dataset.row(4) == Vector<double>{5, 1});
    CHECK(dataset.row(5) == Vector<double>{6, 3});
    CHECK(dataset.row(6) == Vector<double>{7, 6});
    CHECK(dataset.row(7) == Vector<double>{8, 7});
    CHECK(dataset.row(8) == Vector<double>{9, 9});
}

TEST_CASE_METHOD(fixture, "Dataset::size") {
    CHECK(dataset.size() == 9);
}

TEST_CASE_METHOD(fixture, "Dataset::empty") {
    CHECK(!dataset.empty());

    Dataset empty_dataset;
    CHECK(empty_dataset.empty());
}

TEST_CASE("Dataset::save") {
    Dataset dataset({
        std::vector<double>{0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10}, 
        std::vector<double>{1,    2,    3,    4,    5,    6,    7,    8,    9,    10}
    });

    std::string path = "temp/dataset/save.dat";
    dataset.save(path);
    Dataset loaded_dataset(path);
    CHECK(dataset == loaded_dataset);
}

TEST_CASE_METHOD(fixture, "Dataset::col_names") {
    dataset.set_col_names({"a", "b"});
    CHECK(dataset.get_col_names() == std::vector<std::string>{"a", "b"});
}

TEST_CASE("Dataset::interpolate") {    
    SECTION("simple") {
        Dataset data({
            std::vector<double>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, 
            std::vector<double>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}
        });

        data = data.interpolate(1);
        REQUIRE(data.size() == 18);
        CHECK_THAT(data.x(0), Catch::Matchers::WithinAbs(1, 1e-6));
        CHECK_THAT(data.y(0), Catch::Matchers::WithinAbs(1, 1e-6));

        CHECK_THAT(data.x(1), Catch::Matchers::WithinAbs(1.5, 1e-6));
        CHECK_THAT(data.y(1), Catch::Matchers::WithinAbs(1.5, 1e-6));

        CHECK_THAT(data.x(2), Catch::Matchers::WithinAbs(2, 1e-6));
        CHECK_THAT(data.y(2), Catch::Matchers::WithinAbs(2, 1e-6));

        CHECK_THAT(data.x(3), Catch::Matchers::WithinAbs(2.5, 1e-6));
        CHECK_THAT(data.y(3), Catch::Matchers::WithinAbs(2.5, 1e-6));

        CHECK_THAT(data.x(4), Catch::Matchers::WithinAbs(3, 1e-6));
        CHECK_THAT(data.y(4), Catch::Matchers::WithinAbs(3, 1e-6));

        CHECK_THAT(data.x(5), Catch::Matchers::WithinAbs(3.5, 1e-6));
        CHECK_THAT(data.y(5), Catch::Matchers::WithinAbs(3.5, 1e-6));
    }

    SECTION("sine") {
        std::vector<double> x, y;
        for (double xx = 0; xx < 2*M_PI; xx += 0.05) {
            x.push_back(xx);
            y.push_back(sin(xx));
        }

        Dataset data({x, y});
        data = data.interpolate(5);
        for (unsigned int i = 0; i < data.size(); i++) {
            CHECK_THAT(data.y(i), Catch::Matchers::WithinAbs(sin(data.x(i)), 1e-3));
        }
    }

    SECTION("vector interpolation") {
        std::vector<double> x1, y1, x2;
        for (double xx = 0; xx < 2*M_PI; xx += 0.05) {
            x1.push_back(xx);
            y1.push_back(sin(xx));
            x2.push_back(xx + 0.025);
        }

        Dataset data1({x1, y1});
        auto data2 = data1.interpolate(x2);
        for (unsigned int i = 0; i < data2.size(); i++) {
            CHECK_THAT(data2.y(i), Catch::Matchers::WithinAbs(sin(data2.x(i)), 1e-3));
        }
    }

    SECTION("single values") {
        std::vector<double> x1, y1;
        for (double xx = 0; xx < 2*M_PI; xx += 0.05) {
            x1.push_back(xx);
            y1.push_back(sin(xx));
        }

        Dataset data1({x1, y1});
        for (unsigned int i = 0; i < data1.size(); i++) {
            CHECK_THAT(data1.interpolate_y(data1.x(i)+0.025), Catch::Matchers::WithinAbs(sin(data1.x(i)+0.025), 1e-3));
        }
    }
}

TEST_CASE("Dataset::rolling_average") {
    Dataset data({
        std::vector<double>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, 
        std::vector<double>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, 
        std::vector<double>{1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
    });

    SECTION("half_moving_average") {
        SECTION("3") {
            Dataset res = data.rolling_average(3);
            REQUIRE(res.size() == 10);
            CHECK_THAT(res.x(0), Catch::Matchers::WithinAbs(1, 1e-6));
            CHECK_THAT(res.y(0), Catch::Matchers::WithinAbs(1, 1e-6));

            CHECK_THAT(res.x(1), Catch::Matchers::WithinAbs(2, 1e-6));
            CHECK_THAT(res.y(1), Catch::Matchers::WithinAbs((1./2 + 2 + 3./2)/2, 1e-6));

            CHECK_THAT(res.x(2), Catch::Matchers::WithinAbs(3, 1e-6));
            CHECK_THAT(res.y(2), Catch::Matchers::WithinAbs((2./2 + 3 + 4./2)/2, 1e-6));

            CHECK_THAT(res.x(3), Catch::Matchers::WithinAbs(4, 1e-6));
            CHECK_THAT(res.y(3), Catch::Matchers::WithinAbs((3./2 + 4 + 5./2)/2, 1e-6));

            CHECK_THAT(res.x(4), Catch::Matchers::WithinAbs(5, 1e-6));
            CHECK_THAT(res.y(4), Catch::Matchers::WithinAbs((4./2 + 5 + 6./2)/2, 1e-6));

            CHECK_THAT(res.x(5), Catch::Matchers::WithinAbs(6, 1e-6));
            CHECK_THAT(res.y(5), Catch::Matchers::WithinAbs((5./2 + 6 + 7./2)/2, 1e-6));

            CHECK_THAT(res.x(6), Catch::Matchers::WithinAbs(7, 1e-6));
            CHECK_THAT(res.y(6), Catch::Matchers::WithinAbs((6./2 + 7 + 8./2)/2, 1e-6));

            CHECK_THAT(res.x(7), Catch::Matchers::WithinAbs(8, 1e-6));
            CHECK_THAT(res.y(7), Catch::Matchers::WithinAbs((7./2 + 8 + 9./2)/2, 1e-6));

            CHECK_THAT(res.x(8), Catch::Matchers::WithinAbs(9, 1e-6));
            CHECK_THAT(res.y(8), Catch::Matchers::WithinAbs((8./2 + 9 + 10./2)/2, 1e-6));

            CHECK_THAT(res.x(9), Catch::Matchers::WithinAbs(10, 1e-6));
            CHECK_THAT(res.y(9), Catch::Matchers::WithinAbs(10, 1e-6));
        }

        SECTION("5") {
            Dataset res = data.rolling_average(5);
            REQUIRE(res.size() == 10);
            CHECK_THAT(res.x(0), Catch::Matchers::WithinAbs(1, 1e-6));
            CHECK_THAT(res.y(0), Catch::Matchers::WithinAbs(1, 1e-6));

            CHECK_THAT(res.x(1), Catch::Matchers::WithinAbs(2, 1e-6));
            CHECK_THAT(res.y(1), Catch::Matchers::WithinAbs((1./2 + 2 + 3./2)/2, 1e-6));

            CHECK_THAT(res.x(2), Catch::Matchers::WithinAbs(3, 1e-6));
            CHECK_THAT(res.y(2), Catch::Matchers::WithinAbs((1./4 + 2./2 + 3 + 4./2 + 5./4)/2.5, 1e-6));

            CHECK_THAT(res.x(3), Catch::Matchers::WithinAbs(4, 1e-6));
            CHECK_THAT(res.y(3), Catch::Matchers::WithinAbs((2./4 + 3./2 + 4 + 5./2 + 6./4)/2.5, 1e-6));

            CHECK_THAT(res.x(4), Catch::Matchers::WithinAbs(5, 1e-6));
            CHECK_THAT(res.y(4), Catch::Matchers::WithinAbs((3./4 + 4./2 + 5 + 6./2 + 7./4)/2.5, 1e-6));

            CHECK_THAT(res.x(5), Catch::Matchers::WithinAbs(6, 1e-6));
            CHECK_THAT(res.y(5), Catch::Matchers::WithinAbs((4./4 + 5./2 + 6 + 7./2 + 8./4)/2.5, 1e-6));

            CHECK_THAT(res.x(6), Catch::Matchers::WithinAbs(7, 1e-6));
            CHECK_THAT(res.y(6), Catch::Matchers::WithinAbs((5./4 + 6./2 + 7 + 8./2 + 9./4)/2.5, 1e-6));

            CHECK_THAT(res.x(7), Catch::Matchers::WithinAbs(8, 1e-6));
            CHECK_THAT(res.y(7), Catch::Matchers::WithinAbs((6./4 + 7./2 + 8 + 9./2 + 10./4)/2.5, 1e-6));

            CHECK_THAT(res.x(8), Catch::Matchers::WithinAbs(9, 1e-6));
            CHECK_THAT(res.y(8), Catch::Matchers::WithinAbs((8./2 + 9 + 10./2)/2, 1e-6));

            CHECK_THAT(res.x(9), Catch::Matchers::WithinAbs(10, 1e-6));
            CHECK_THAT(res.y(9), Catch::Matchers::WithinAbs(10, 1e-6));
        }
    }
}

TEST_CASE_METHOD(fixture, "Dataset::append") {
    auto d2 = dataset;
    dataset.append(d2);
    REQUIRE(dataset.size() == 18);
    for (unsigned int i = 0; i < 9; ++i) {
        CHECK(dataset.x(i) == dataset.x(i+9));
        CHECK(dataset.y(i) == dataset.y(i+9));
    }

    d2 = dataset;
    dataset.append(d2);
    REQUIRE(dataset.size() == 36);
    for (unsigned int i = 0; i < 18; ++i) {
        CHECK(dataset.x(i) == dataset.x(i+18));
        CHECK(dataset.y(i) == dataset.y(i+18));
    }
}

TEST_CASE_METHOD(fixture, "Dataset::limit_x") {
    SECTION("simple") {
        Limit limit(0.5, 5);
        dataset.limit_x(limit);
        for (unsigned int i = 0; i < dataset.size(); i++) {
            CHECK(limit.min <= dataset.x(i));
            CHECK(dataset.x(i) <= limit.max);
        }
    }

    SECTION("real data") {
        Dataset data("test/files/2epe.dat");

        unsigned int start = 0;
        while (data.x(start) < 0.01) {            
            start++;
        }

        unsigned int end = data.size()-1;
        while (0.3 < data.x(end)) {
            end--;
        }

        auto data_limited = data;
        data_limited.limit_x(0.01, 0.3);
        REQUIRE(data_limited.size() == end-start+1);
        for (unsigned int i = 0; i < data_limited.size(); i++) {
            CHECK(data_limited.x(i) == data.x(i+start));
            CHECK(data_limited.y(i) == data.y(i+start));
        }
    }
}

TEST_CASE_METHOD(fixture, "Dataset::limit_y") {
    Limit limit(0.5, 5);
    dataset.limit_y(limit);
    for (unsigned int i = 0; i < dataset.size(); i++) {
        CHECK(limit.min <= dataset.y(i));
        CHECK(dataset.y(i) <= limit.max);
    }
}

TEST_CASE_METHOD(fixture, "Dataset::x") {
    CHECK(dataset.x() == dataset.col(0));
}

TEST_CASE_METHOD(fixture, "Dataset::y") {
    CHECK(dataset.y() == dataset.col(1));
}

TEST_CASE("Dataset::sort_x") {
    std::vector<double> x = {1, 0, 4, 3, 2};
    std::vector<double> y = {10, 0, 40, 30, 20};
    Dataset data({x, y});

    data.sort_x();
    CHECK(data.x() == std::vector<double>({0, 1, 2, 3, 4}));
    CHECK(data.y() == std::vector<double>({0, 10, 20, 30, 40}));
}

std::string generate_SASDJG5_dataset();
TEST_CASE("Dataset::find_minima") {
    // SECTION("empty") {
    //     Dataset data;
    //     std::vector<unsigned int> minima = data.find_minima();
    //     REQUIRE(minima.empty());
    // }

    // SECTION("simple") {
    //     SECTION("single") {
    //         std::vector<double> x = {1,  2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
    //         std::vector<double> y = {10, 9, 8, 7, 6, 5, 4, 3, 2, 1,  2,  3,  4,  5,  6,  7,  8,  9,  10};
    //         Dataset data({x, y});
    //         std::vector<unsigned int> minima = data.find_minima();
    //         REQUIRE(minima == std::vector<unsigned int>{9});
    //     }

        // SECTION("multiple") {
        //     std::vector<double> x = {1,  2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
        //     std::vector<double> y = {10, 9, 8, 7, 6, 7, 8, 9, 8, 7,  6,  5,  4,  5,  6,  7,  8,  9,  10};
        //     Dataset data({x, y});
        //     std::vector<unsigned int> minima = data.find_minima();
        //     REQUIRE(minima == std::vector<unsigned int>{4, 12});
        // }

    //     SECTION("at endpoints") {
    //         std::vector<double> x = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    //         std::vector<double> y = {1, 2, 3, 4, 5, 4, 3, 2, 1};
    //         Dataset data({x, y});
    //         std::vector<unsigned int> minima = data.find_minima();
    //         REQUIRE(minima == std::vector<unsigned int>{0, 8});
    //     }
    // }

    // SECTION("sinusoidal noise") {
    //     for (unsigned int i = 0; i < 10; i++) {
    //         std::vector<double> x;
    //         std::vector<double> y;
    //         for (double xx = -10, dx = 0.1; xx <= 10; xx += dx) {
    //             x.push_back(xx);
    //             y.push_back(0.5*xx*xx + 2*std::sin(2*xx)*std::rand()/RAND_MAX);
    //         }
    //         Dataset data({x, y});
    //         data = data.rolling_average(7);
    //         std::vector<unsigned int> minima = data.find_minima(1, 0.05);
    //         REQUIRE(minima.size() < 5);
    //         for (unsigned int i = 0; i < minima.size()-1; i++) {
    //             CHECK(-2 < x[minima[i]]);
    //             CHECK(x[minima[i]] < 2);
    //         }
    //     }
    // }

    SECTION("actual data") {
        auto file = generate_SASDJG5_dataset();
        Dataset data(file);
        data = data.rolling_average(5);

        // find all minima
        auto minima = data.find_minima();
        REQUIRE(minima.size() == 5);
        CHECK(minima[0] == 4);
        CHECK(minima[1] == 41);
        CHECK(minima[2] == 48);
        CHECK(minima[3] == 66);
        CHECK(minima[4] == 73);

        // remove one by increasing the prominence
        minima = data.find_minima(1, 0.05);
        REQUIRE(minima.size() == 4);
        CHECK(minima[0] == 4);
        CHECK(minima[1] == 41);
        CHECK(minima[2] == 48);
        CHECK(minima[3] == 66);

        // remove another by further increasing prominence
        minima = data.find_minima(1, 0.15);
        REQUIRE(minima.size() == 3);
        CHECK(minima[0] == 4);
        CHECK(minima[1] == 41);
        CHECK(minima[2] == 66);

        // remove all but 2 by using a large prominence
        minima = data.find_minima(1, 0.8);
        REQUIRE(minima.size() == 2);
        CHECK(minima[0] == 4);
        CHECK(minima[1] == 41);

    //     SECTION("find_maxima") {
    //         auto maxima = data.find_maxima();
    //         REQUIRE(2 <= maxima.size());
    //         CHECK(maxima[0] == 9);
    //         CHECK(maxima[1] == 19);
    //     }
    }
}

std::string generate_SASDJG5_dataset() {
    io::File file("temp/test/dataset/SASDJG5.dat");
    if (file.exists()) {return file;}
    file.create();

    std::string data = 
    "1.54196666e-03   7.54563335e+03   0.00000000e+00\n" 
    "1.79474393e-03   7.36448607e+03   0.00000000e+00\n" 
    "2.04752120e-03   7.26108566e+03   0.00000000e+00\n" 
    "2.30029847e-03   7.20216108e+03   0.00000000e+00\n" 
    "2.55307574e-03   7.18847353e+03   0.00000000e+00\n" 
    "2.80585301e-03   7.18884069e+03   0.00000000e+00\n" 
    "3.05863028e-03   7.22114198e+03   0.00000000e+00\n" 
    "3.31140755e-03   7.25230012e+03   0.00000000e+00\n" 
    "3.56418482e-03   7.32054039e+03   0.00000000e+00\n" 
    "3.81696209e-03   7.39549653e+03   0.00000000e+00\n" 
    "4.06973936e-03   7.55560453e+03   0.00000000e+00\n" 
    "4.32251663e-03   7.71682064e+03   0.00000000e+00\n" 
    "4.57529391e-03   7.94209171e+03   0.00000000e+00\n" 
    "4.82807118e-03   8.23176279e+03   0.00000000e+00\n" 
    "5.08084845e-03   8.56771713e+03   0.00000000e+00\n" 
    "5.33362572e-03   8.87773916e+03   0.00000000e+00\n" 
    "5.58640299e-03   9.28444071e+03   0.00000000e+00\n" 
    "5.83918026e-03   9.66772687e+03   0.00000000e+00\n" 
    "6.09195753e-03   1.01475205e+04   0.00000000e+00\n" 
    "6.34473480e-03   1.06754464e+04   0.00000000e+00\n" 
    "6.59751207e-03   1.11850792e+04   0.00000000e+00\n" 
    "6.85028934e-03   1.16868057e+04   0.00000000e+00\n" 
    "7.10306661e-03   1.22415133e+04   0.00000000e+00\n" 
    "7.35584389e-03   1.28070981e+04   0.00000000e+00\n" 
    "7.60862116e-03   1.34043749e+04   0.00000000e+00\n" 
    "7.86139843e-03   1.35158485e+04   0.00000000e+00\n" 
    "8.11417570e-03   1.39970611e+04   0.00000000e+00\n" 
    "8.36695297e-03   1.40636833e+04   0.00000000e+00\n" 
    "8.61973024e-03   1.41430466e+04   0.00000000e+00\n" 
    "8.87250751e-03   1.37648877e+04   0.00000000e+00\n" 
    "9.12528478e-03   1.32558786e+04   0.00000000e+00\n" 
    "9.37806205e-03   1.26512085e+04   0.00000000e+00\n" 
    "9.63083932e-03   1.19439015e+04   0.00000000e+00\n" 
    "9.88361659e-03   1.12143305e+04   0.00000000e+00\n" 
    "1.01363939e-02   1.02245963e+04   0.00000000e+00\n" 
    "1.03891711e-02   9.41102157e+03   0.00000000e+00\n" 
    "1.06419484e-02   8.07823600e+03   0.00000000e+00\n" 
    "1.08947257e-02   7.52041074e+03   0.00000000e+00\n" 
    "1.11475029e-02   6.34829824e+03   0.00000000e+00\n" 
    "1.14002802e-02   5.44554795e+03   0.00000000e+00\n" 
    "1.16530575e-02   5.46443444e+03   0.00000000e+00\n" 
    "1.19058348e-02   5.28319823e+03   0.00000000e+00\n" 
    "1.21586120e-02   4.93155586e+03   0.00000000e+00\n" 
    "1.24113893e-02   5.72876316e+03   0.00000000e+00\n" 
    "1.26641666e-02   6.21435772e+03   0.00000000e+00\n" 
    "1.29169438e-02   6.81131801e+03   0.00000000e+00\n" 
    "1.31697211e-02   6.16321675e+03   0.00000000e+00\n" 
    "1.34224984e-02   5.49700207e+03   0.00000000e+00\n" 
    "1.36752757e-02   6.23447959e+03   0.00000000e+00\n" 
    "1.39280529e-02   5.81354805e+03   0.00000000e+00\n" 
    "1.41808302e-02   6.09721155e+03   0.00000000e+00\n" 
    "1.44336075e-02   7.43127854e+03   0.00000000e+00\n" 
    "1.46863847e-02   7.82821797e+03   0.00000000e+00\n" 
    "1.49391620e-02   9.22595753e+03   0.00000000e+00\n" 
    "1.51919393e-02   7.50807824e+03   0.00000000e+00\n" 
    "1.54447166e-02   9.95165150e+03   0.00000000e+00\n" 
    "1.56974938e-02   1.06986044e+04   0.00000000e+00\n" 
    "1.59502711e-02   1.05103991e+04   0.00000000e+00\n" 
    "1.62030484e-02   1.18825706e+04   0.00000000e+00\n" 
    "1.64558256e-02   1.31441838e+04   0.00000000e+00\n" 
    "1.67086029e-02   1.40661017e+04   0.00000000e+00\n" 
    "1.69613802e-02   1.21325725e+04   0.00000000e+00\n" 
    "1.72141574e-02   1.79476353e+04   0.00000000e+00\n" 
    "1.74669347e-02   1.69114349e+04   0.00000000e+00\n" 
    "1.77197120e-02   1.19884817e+04   0.00000000e+00\n" 
    "1.79724893e-02   1.25151799e+04   0.00000000e+00\n" 
    "1.82252665e-02   1.53953789e+04   0.00000000e+00\n" 
    "1.84780438e-02   1.61259638e+04   0.00000000e+00\n" 
    "1.87308211e-02   1.26381116e+04   0.00000000e+00\n" 
    "1.89835983e-02   1.61748750e+04   0.00000000e+00\n" 
    "1.92363756e-02   2.04171060e+04   0.00000000e+00\n" 
    "1.94891529e-02   1.94805321e+04   0.00000000e+00\n" 
    "1.97419302e-02   1.78842601e+04   0.00000000e+00\n" 
    "1.99947074e-02   2.07731569e+04   0.00000000e+00\n" 
    "2.02474847e-02   2.11909101e+04   0.00000000e+00\n" 
    "2.05002620e-02   1.67068997e+04   0.00000000e+00\n" 
    "2.07530392e-02   2.04606532e+04   0.00000000e+00\n" 
    "2.12585938e-02   2.16886419e+04   0.00000000e+00\n" 
    "2.17641483e-02   2.43663580e+04   0.00000000e+00\n" 
    "2.20169256e-02   2.36571727e+04   0.00000000e+00\n" 
    "2.32808119e-02   2.41054892e+04   0.00000000e+00\n";

    std::ofstream out(file);
    out << data;
    out.close();
    return file;
}