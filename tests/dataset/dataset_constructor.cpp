#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <dataset/detail/DATReader.h>
#include <dataset/detail/XVGReader.h>
#include <dataset/Dataset.h>
#include <math/Vector.h>
#include <settings/GeneralSettings.h>
#include <settings/HistogramSettings.h>

#include <fstream>

using namespace ausaxs;

TEST_CASE("DATReader::construct") {
    settings::general::verbose = false;
    io::File test_file = "temp/tests/dataset/dat_test.dat";
    test_file.directory().create();

    SECTION("simple contents") {
        std::string test_file_contents = 
            "x y z\n"
            "0.1 1 10\n"
            "0.2 2 20\n"
            "0.3 3 30\n";

        std::ofstream file(test_file);
        file << test_file_contents;
        file.close();

        // default M
        auto data = detail::DATReader().construct(test_file, 0);
        REQUIRE(data->size_cols() == 3);
        CHECK(data->col(0) == Vector<double>({0.1, 0.2, 0.3}));
        CHECK(data->col(1) == Vector<double>({1, 2, 3}));
        CHECK(data->col(2) == Vector<double>({10, 20, 30}));

        // M = 1
        auto data1 = detail::DATReader().construct(test_file, 1);
        REQUIRE(data1->size_cols() == 1);
        CHECK(data1->col(0) == Vector<double>({0.1, 0.2, 0.3}));

        // M = 2
        auto data2 = detail::DATReader().construct(test_file, 2);
        REQUIRE(data2->size_cols() == 2);
        CHECK(data2->col(0) == Vector<double>({0.1, 0.2, 0.3}));
        CHECK(data2->col(1) == Vector<double>({1, 2, 3}));

        // M = 3
        auto data3 = detail::DATReader().construct(test_file, 3);
        REQUIRE(data3->size_cols() == 3);
        CHECK(data3->col(0) == Vector<double>({0.1, 0.2, 0.3}));
        CHECK(data3->col(1) == Vector<double>({1, 2, 3}));
        CHECK(data3->col(2) == Vector<double>({10, 20, 30}));
    }

    SECTION("weird contents") {
        std::string test_file_contents = 
            "x y z\n"
            "0.1 1 10 100\n"
            "0.11 1.1 11\n"
            "0.12 1.2\n"
            "skip me\n"
            "0.2 2 20 200\n"
            "0.3 3 30 300\n"
            "0.4 4 40 400\n";

        {
            std::ofstream file(test_file);
            file << test_file_contents;
            file.close();

            // default M
            auto data = detail::DATReader().construct(test_file, 0);
            REQUIRE(data->size_cols() == 4);
            CHECK(data->col(0) == Vector<double>({0.1, 0.2, 0.3, 0.4}));
            CHECK(data->col(1) == Vector<double>({1,   2,   3,   4}));
            CHECK(data->col(2) == Vector<double>({10,  20,  30,  40}));
            CHECK(data->col(3) == Vector<double>({100, 200, 300, 400}));

            // M = 1
            auto data1 = detail::DATReader().construct(test_file, 1);
            REQUIRE(data1->size_cols() == 1);
            CHECK(data1->col(0) == Vector<double>({0.1, 0.2, 0.3, 0.4}));

            // M = 2
            auto data2 = detail::DATReader().construct(test_file, 2);
            REQUIRE(data2->size_cols() == 2);
            CHECK(data2->col(0) == Vector<double>({0.1, 0.2, 0.3, 0.4}));
            CHECK(data2->col(1) == Vector<double>({1,   2,   3,   4}));

            // M = 3
            auto data3 = detail::DATReader().construct(test_file, 3);
            REQUIRE(data3->size_cols() == 3);
            CHECK(data3->col(0) == Vector<double>({0.1, 0.2, 0.3, 0.4}));
            CHECK(data3->col(1) == Vector<double>({1,   2,   3,   4}));
            CHECK(data3->col(2) == Vector<double>({10,  20,  30,  40}));

            // M = 4
            auto data4 = detail::DATReader().construct(test_file, 4);
            REQUIRE(data4->size_cols() == 4);
            CHECK(data4->col(0) == Vector<double>({0.1, 0.2, 0.3, 0.4}));
            CHECK(data4->col(1) == Vector<double>({1,   2,   3,   4}));
            CHECK(data4->col(2) == Vector<double>({10,  20,  30,  40}));
            CHECK(data4->col(3) == Vector<double>({100, 200, 300, 400}));
        }
    }
}

TEST_CASE("Dataset: different unit") {
    settings::general::verbose = false;

    Dataset dataset({
        std::vector<double>{0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10}, 
        std::vector<double>{1,    2,    3,    4,    5,    6,    7,    8,    9,    10}
    });

    SECTION("specified in file") {
        settings::axes::qmax = 1;

        dataset.save(          "temp/tests/dataset/save.dat", "[nm]");
        Dataset loaded_dataset("temp/tests/dataset/save.dat");
        REQUIRE(loaded_dataset.size() == dataset.size());
        for (unsigned int i = 0; i < dataset.size(); i++) {
            REQUIRE_THAT(loaded_dataset.x(i)*10, Catch::Matchers::WithinAbs(dataset.x(i), 1e-6));
            REQUIRE(loaded_dataset.y(i) == dataset.y(i));
        }
    }

    SECTION("specified by setting") {
        settings::general::input_q_unit = settings::general::QUnit::NM;
        dataset.save(          "temp/tests/dataset/save.dat");
        Dataset loaded_dataset("temp/tests/dataset/save.dat");
        REQUIRE(loaded_dataset.size() == dataset.size());
        for (unsigned int i = 0; i < dataset.size(); i++) {
            REQUIRE_THAT(loaded_dataset.x(i)*10, Catch::Matchers::WithinAbs(dataset.x(i), 1e-6));
            REQUIRE(loaded_dataset.y(i) == dataset.y(i));
        }
        settings::general::input_q_unit = settings::general::QUnit::A;
    }
}

auto vec_approx = [](const auto& v1, const auto& v2) {
    REQUIRE(v1.size() == v2.size());
    for (unsigned int i = 0; i < v1.size(); i++) {
        CHECK_THAT(v1[i], Catch::Matchers::WithinAbs(v2[i], 1e-6));
    }
};

TEST_CASE("XVGReader::construct") {
    settings::general::verbose = false;
    io::File test_file = "temp/tests/dataset/dat_test.dat";
    test_file.directory().create();

    SECTION("simple contents") {
        std::string test_file_contents = 
            "x y z\n"
            "0.1 1 10\n"
            "0.2 2 20\n"
            "0.3 3 30\n";

        {
            std::ofstream file(test_file);
            file << test_file_contents;
            file.close();

            // default M
            auto data = detail::XVGReader().construct(test_file, 0);
            REQUIRE(data->size_cols() == 3);
            vec_approx(data->col(0), std::vector<double>({0.01, 0.02, 0.03}));
            vec_approx(data->col(1), std::vector<double>({1, 2, 3}));
            vec_approx(data->col(2), std::vector<double>({10, 20, 30}));

            // M = 1
            auto data1 = detail::XVGReader().construct(test_file, 1);
            REQUIRE(data1->size_cols() == 1);
            vec_approx(data1->col(0), std::vector<double>({0.01, 0.02, 0.03}));

            // M = 2
            auto data2 = detail::XVGReader().construct(test_file, 2);
            REQUIRE(data2->size_cols() == 2);
            vec_approx(data2->col(0), std::vector<double>({0.01, 0.02, 0.03}));
            vec_approx(data2->col(1), std::vector<double>({1,    2,    3}));

            // M = 3
            auto data3 = detail::XVGReader().construct(test_file, 3);
            REQUIRE(data3->size_cols() == 3);
            vec_approx(data3->col(0), std::vector<double>({0.01, 0.02, 0.03}));
            vec_approx(data3->col(1), std::vector<double>({1,    2,    3}));
            vec_approx(data3->col(2), std::vector<double>({10,   20,   30}));
        }
    }

    SECTION("weird contents") {
        std::string test_file_contents = 
            "x y z\n"
            "0.1 1 10 100\n"
            "0.11 1.1 11\n"
            "0.12 1.2\n"
            "skip me\n"
            "0.2 2 20 200\n"
            "0.3 3 30 300\n"
            "0.4 4 40 400\n";

        {
            std::ofstream file(test_file);
            file << test_file_contents;
            file.close();

            // default M
            auto data = detail::XVGReader().construct(test_file, 0);
            REQUIRE(data->size_cols() == 4);
            vec_approx(data->col(0), std::vector<double>({0.01, 0.02, 0.03, 0.04}));
            vec_approx(data->col(1), std::vector<double>({1,    2,    3,    4}));
            vec_approx(data->col(2), std::vector<double>({10,   20,   30,   40}));
            vec_approx(data->col(3), std::vector<double>({100,  200,  300,  400}));

            // M = 1
            auto data1 = detail::XVGReader().construct(test_file, 1);
            REQUIRE(data1->size_cols() == 1);
            vec_approx(data1->col(0), std::vector<double>({0.01, 0.02, 0.03, 0.04}));

            // M = 2
            auto data2 = detail::XVGReader().construct(test_file, 2);
            REQUIRE(data2->size_cols() == 2);
            vec_approx(data2->col(0), std::vector<double>({0.01, 0.02, 0.03, 0.04}));
            vec_approx(data2->col(1), std::vector<double>({1,    2,    3,    4}));

            // M = 3
            auto data3 = detail::XVGReader().construct(test_file, 3);
            REQUIRE(data3->size_cols() == 3);
            vec_approx(data3->col(0), std::vector<double>({0.01, 0.02, 0.03, 0.04}));
            vec_approx(data3->col(1), std::vector<double>({1,    2,    3,    4}));
            vec_approx(data3->col(2), std::vector<double>({10,   20,   30,   40}));

            // M = 4
            auto data4 = detail::XVGReader().construct(test_file, 4);
            REQUIRE(data4->size_cols() == 4);
            vec_approx(data4->col(0), std::vector<double>({0.01, 0.02, 0.03, 0.04}));
            vec_approx(data4->col(1), std::vector<double>({1,    2,    3,    4}));
            vec_approx(data4->col(2), std::vector<double>({10,   20,   30,   40}));
            vec_approx(data4->col(3), std::vector<double>({100,  200,  300,  400}));
        }
    }
}