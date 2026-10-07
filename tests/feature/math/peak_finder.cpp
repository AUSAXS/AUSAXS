#include <catch2/catch_test_macros.hpp>

#include <dataset/Dataset.h>
#include <math/PeakFinder.h>
#include <settings/GeneralSettings.h>

#include <support/temp_file.h>

#include <cmath>
#include <cstdlib>

using namespace ausaxs;

static test::TempFile generate_SASDJG5_dataset();
TEST_CASE("math::find_minima") {
    settings::general::verbose = false;
    SECTION("simple") {
        SECTION("single") {
            std::vector<double> x = {1,  2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
            std::vector<double> y = {10, 9, 8, 7, 6, 5, 4, 3, 2, 1,  2,  3,  4,  5,  6,  7,  8,  9,  10};
            Dataset data({x, y});
            std::vector<int> minima = math::find_minima(data.x(), data.y(), 0, 0);
            REQUIRE(minima == std::vector<int>{9});
        }

        SECTION("multiple") {
            std::vector<double> x = {1,  2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
            std::vector<double> y = {10, 9, 8, 7, 6, 7, 8, 9, 8, 7,  6,  5,  4,  5,  6,  7,  8,  9,  10};
            Dataset data({x, y});
            std::vector<int> minima = math::find_minima(data.x(), data.y(), 0, 0);
            REQUIRE(minima == std::vector<int>{4, 12});
        }

        SECTION("at endpoints") {
            std::vector<double> x = {1, 2, 3, 4, 5, 6, 7, 8, 9};
            std::vector<double> y = {1, 2, 3, 4, 5, 4, 3, 2, 1};
            Dataset data({x, y});
            std::vector<int> minima = math::find_minima(data.x(), data.y(), 0, 0);
            REQUIRE(minima == std::vector<int>{0, 8});
        }
    }

    SECTION("sinusoidal noise") {
        for (int i = 0; i < 10; i++) {
            std::vector<double> x;
            std::vector<double> y;
            for (double xx = -10, dx = 0.1; xx <= 10; xx += dx) {
                x.push_back(xx);
                y.push_back(0.5*xx*xx + 2*std::sin(2*xx)*std::rand()/RAND_MAX);
            }
            Dataset data({x, y});

            data = data.rolling_average(7);
            std::vector<int> minima = math::find_minima(data.x(), data.y(), 1, 0.05);
            REQUIRE(minima.size() < 5);
            for (int j = 0; j < static_cast<int>(minima.size())-1; j++) {
                CHECK(-2 < x[minima[j]]);
                CHECK(x[minima[j]] < 2);
            }
        }
    }

    SECTION("actual data") {
        auto file = generate_SASDJG5_dataset();
        Dataset data(file);
        data = data.rolling_average(5);

        // find all minima
        auto minima = math::find_minima(data.x(), data.y(), 0, 0);
        REQUIRE(minima.size() == 5);
        CHECK(minima[0] == 4);
        CHECK(minima[1] == 41);
        CHECK(minima[2] == 48);
        CHECK(minima[3] == 66);
        CHECK(minima[4] == 73);

        // remove one by increasing the prominence
        minima = math::find_minima(data.x(), data.y(), 1, 0.05);
        REQUIRE(minima.size() == 4);
        CHECK(minima[0] == 4);
        CHECK(minima[1] == 41);
        CHECK(minima[2] == 48);
        CHECK(minima[3] == 66);

        // remove another by further increasing prominence
        minima = math::find_minima(data.x(), data.y(), 1, 0.15);
        REQUIRE(minima.size() == 3);
        CHECK(minima[0] == 4);
        CHECK(minima[1] == 41);
        CHECK(minima[2] == 66);

        // remove all but 2 by using a large prominence
        minima = math::find_minima(data.x(), data.y(), 1, 0.8);
        REQUIRE(minima.size() == 2);
        CHECK(minima[0] == 4);
        CHECK(minima[1] == 41);
    }
}

TEST_CASE("math::find_minima: maxima of inverted data") {
    settings::general::verbose = false;
    SECTION("actual data") {
        auto file = generate_SASDJG5_dataset();
        Dataset data(file);
        data = data.rolling_average(5);

        // find all maxima
        auto minima = math::find_minima(data.x(), -data.y(), 0, 0);
        REQUIRE(minima.size() == 5);
        CHECK(minima[0] == 0);
        CHECK(minima[1] == 27);
        CHECK(minima[2] == 45);
        CHECK(minima[3] == 62);
        CHECK(minima[4] == 80);

        // remove one by increasing the prominence
        minima = math::find_minima(data.x(), -data.y(), 1, 0.05);
        REQUIRE(minima.size() == 3);
        CHECK(minima[0] == 27);
        CHECK(minima[1] == 62);
        CHECK(minima[2] == 80);

        // remove all but 2 by using a large prominence
        minima = math::find_minima(data.x(), -data.y(), 1, 0.8);
        REQUIRE(minima.size() == 2);
        CHECK(minima[0] == 27);
        CHECK(minima[1] == 80);
    }
}

test::TempFile generate_SASDJG5_dataset() {
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

    return test::TempFile(".dat", data);
}