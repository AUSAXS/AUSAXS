#include "settings/GeneralSettings.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <dataset/SimpleDataset.h>
#include <math/Statistics.h>

using namespace ausaxs;

struct fixture {
    SimpleDataset dataset = {{1, 2, 3}, {4, 5, 6}};
};

TEST_CASE("SimpleDataset::SimpleDataset") {
    settings::general::verbose = false;
    SECTION("io::ExistingFile") {
        SimpleDataset dataset("tests/files/2epe.dat");
        CHECK(dataset.size() == 104);
        CHECK(dataset.size_rows() == 104);
        CHECK(dataset.size_cols() == 3);
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
        std::string path = "temp/dataset/save.dat";
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
        SimpleDataset data("tests/files/2epe.dat");
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
    settings::general::verbose = false;
    SECTION("with real data") {
        SimpleDataset data("tests/files/2epe.dat");
        CHECK(data.size() == 104);
        data.reduce(25, true);
        CHECK(data.size() <= 25);
    }

    SECTION("explicit") {
        std::vector<double> x = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        std::vector<double> y = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        SimpleDataset data(x, y);

        data.reduce(6, true);
        CHECK(data.size_rows() == 6);
        CHECK(data.size_cols() == 3);
        CHECK(data.x() == std::vector{1, 2, 3, 5, 7, 10});
    }

    SECTION("10000 to 100") {
        std::vector<double> x(10000);
        std::vector<double> y(10000);
        for (unsigned int i = 0; i < x.size(); i++) {
            x[i] = std::pow(10, i*6.0/10000);
            y[i] = i;
        }
        SimpleDataset data(x, y);

        data.reduce(100, true);
        CHECK(data.size_rows() == 100);
        CHECK(data.size_cols() == 3);
        CHECK(data.x() == std::vector{1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6});
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
    SimpleDataset dataset("tests/files/2epe.dat");
    CHECK(dataset.size() == 104);
    dataset.rebin();
    CHECK(dataset.size() < 104);
}
