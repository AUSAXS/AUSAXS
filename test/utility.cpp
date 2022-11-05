#include <catch2/catch_test_macros.hpp>
#include <fitter/FitReporter.h>
#include <plots/PlotOptions.h>
#include <plots/PlotDataset.h>
#include <utility/Utility.h>
#include <plots/all.h>

using std::vector, std::string;

TEST_CASE("plots", "[utility],[manual]") {
    std::vector<double> x = {-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5};
    std::vector<double> y = {-5, -3, 4, 8, 12, 7, 3, 1, -3, -5, -9};
    std::vector<double> yerr = vector<double>(y.size(), 0.5);
    SimpleDataset data(x, y, yerr);

    // SECTION("markers") {
    //     plots::PlotDataset::quick_plot(data, "figures/test/utility/plots/default.png");

    //     data.add_plot_options("markers");
    //     plots::PlotDataset::quick_plot(data, "figures/test/utility/plots/markers.png");

    //     data.add_plot_options("errors");
    //     plots::PlotDataset::quick_plot(data, "figures/test/utility/plots/errors.png");
    // }

    // SECTION("ylimits") {
    //     data.add_plot_options("markers", {{"ylimits", Limit(3, 5)}});
    //     plots::PlotDataset::quick_plot(data, "figures/test/utility/plots/limits_y.png");
    // }

    // SECTION("xlimits") {
    //     data.add_plot_options("markers", {{"xlimits", Limit(0, 5)}});
    //     plots::PlotDataset::quick_plot(data, "figures/test/utility/plots/limits_x.png");
    // }

    // SECTION("log") {
    //     data.limit_y(1, inf);
    //     data.add_plot_options({{"logy", true}});
    //     plots::PlotDataset::quick_plot(data, "figures/test/utility/plots/log.png");
    // }

    SECTION("intensityfit") {
        x =                      {0.001, 0.002, 0.003, 0.009, 0.02, 0.06, 0.1, 0.2, 0.4, 0.7};
        std::vector<double> y1 = {10,    10,    9,     9.2,   8.4,  7,    5,   3,   2,   3};
        std::vector<double> y2 = {11,    10.5,  10,    9.5,   8.5,  7.4,  5.6, 3.3, 2.2, 2};
        std::vector<double> y2err = {0.5, 0.5,  0.5,   0.5,   0.5,  0.5,  0.5, 0.5, 0.5, 0.5};

        std::shared_ptr<Fit> fit = std::make_shared<Fit>();
        fit->figures.data = SimpleDataset(x, y1, y2err);
        fit->figures.intensity = SimpleDataset(x, y2);
        fit->figures.intensity_interpolated = SimpleDataset(x, y2);

        plots::PlotIntensityFit::quick_plot(fit, "figures/test/utility/plots/intensityfit.png");
    }
}

TEST_CASE("fitreporter", "[utility],[manual]") {
    Fit fit;
    fit.status = true;
    fit.fevals = 100;
    fit.fval = 1000;
    fit.dof = 3;
    fit.parameters = {{"a", 2, 0.5}, {"b", 3, 0.4}, {"c", 4, 0.3}};

    SECTION("Single") {
        SECTION("Terminal printing") {
            FitReporter::report(fit);
        }
        SECTION("File printing") {
            FitReporter::save("temp/fitreport1.txt", fit);
        }
    }

    SECTION("Multi") {
        Fit fit2;
        fit.status = false;
        fit2.fevals = 20;
        fit2.fval = 200;
        fit2.dof = 6;
        fit2.parameters = {{"a", 1, 0.1}, {"b", 0.5, 0.3}, {"c", 8, 0.6}};

        vector<Fit> fits = {fit, fit2};
        vector<string> titles = {"First fit", "Second fit"};

        SECTION("Terminal printing") {
            FitReporter::report(fits, titles);   
        }
        SECTION("File printing") {
            FitReporter::save("temp/fitreport2.txt", fits, titles);
        }
    }
}

TEST_CASE("split", "[utility]") {
    SECTION("single delimiter") {
        SECTION("empty") {
            vector<string> result = utility::split("", ',');
            REQUIRE(result.size() == 0);
        }

        SECTION("single") {
            vector<string> result = utility::split("a", ',');
            REQUIRE(result.size() == 1);
            REQUIRE(result[0] == "a");
        }

        SECTION("multiple") {
            vector<string> result = utility::split("a,b,c", ',');
            REQUIRE(result.size() == 3);
            REQUIRE(result[0] == "a");
            REQUIRE(result[1] == "b");
            REQUIRE(result[2] == "c");
        }
    }

    SECTION("multiple delimiters") {
        SECTION("empty") {
            vector<string> result = utility::split("", ",");
            REQUIRE(result.size() == 0);
        }

        SECTION("single") {
            vector<string> result = utility::split("a", ",");
            REQUIRE(result.size() == 1);
            REQUIRE(result[0] == "a");
        }

        SECTION("multiple") {
            vector<string> result = utility::split("a,b,c", ",");
            REQUIRE(result.size() == 3);
            REQUIRE(result[0] == "a");
            REQUIRE(result[1] == "b");
            REQUIRE(result[2] == "c");
        }

        SECTION("multiple with spaces") {
            vector<string> result = utility::split("a, b, c", ", ");
            REQUIRE(result.size() == 3);
            REQUIRE(result[0] == "a");
            REQUIRE(result[1] == "b");
            REQUIRE(result[2] == "c");
        }

        SECTION("multiple with lots of stuff") {
            vector<string> result = utility::split("  ,aaba  ,, ,,bda,  ced,,,  ,", ", ");
            REQUIRE(result.size() == 3);
            REQUIRE(result[0] == "aaba");
            REQUIRE(result[1] == "bda");
            REQUIRE(result[2] == "ced");
        }

        SECTION("actual example") {
            vector<string> result = utility::split("0.009813      	0.00667934    	0.00133647    	1\n\r", ", \t\n\r");
            REQUIRE(result.size() == 4);
            REQUIRE(result[0] == "0.009813");
            REQUIRE(result[1] == "0.00667934");
            REQUIRE(result[2] == "0.00133647");
            REQUIRE(result[3] == "1");
        }
    }
}

TEST_CASE("remove_all", "[utility]") {
    SECTION("empty") {
        REQUIRE(utility::remove_all("", "a") == "");
    }

    SECTION("single") {
        REQUIRE(utility::remove_all("a", "a") == "");
    }

    SECTION("multiple") {
        REQUIRE(utility::remove_all("aaa", "a") == "");
    }

    SECTION("multiple with other") {
        REQUIRE(utility::remove_all("ababab", "a") == "bbb");
    }

    SECTION("multiple with other 2") {
        REQUIRE(utility::remove_all("ababab", "b") == "aaa");
    }

    SECTION("real issues") {
        REQUIRE(utility::remove_all("TITLE \r", " \r") == "TITLE");
    }
}

TEST_CASE("join", "[utility]") {
    SECTION("empty") {
        string result = utility::join({}, ",");
        REQUIRE(result == "");
    }

    SECTION("single") {
        string result = utility::join({"a"}, ",");
        REQUIRE(result == "a");
    }

    SECTION("multiple") {
        string result = utility::join({"a", "b", "c"}, ",");
        REQUIRE(result == "a,b,c");
    }

    SECTION("multiple with spaces") {
        string result = utility::join({"a", " b", " c"}, ",");
        REQUIRE(result == "a, b, c");
    }
}

TEST_CASE("plotoptions", "[utility]") {
    plots::PlotOptions options;

    SECTION("set") {
        options.set({{"xlabel", "label for x"}, {"ylabel", string("label for y")}, {"line", true}, {"errors", true}, {"color", style::color::blue}, {"alpha", 0.3}});
        CHECK(options.xlabel == "label for x");
        CHECK(options.ylabel == "label for y");
        CHECK(options.draw_line == true);
        CHECK(options.draw_errors == true);
        CHECK(options.color == style::color::blue);
        CHECK(options.alpha == 0.3);
    }
}

TEST_CASE("utility", "[utility]") {
    string s = "result should be 92 alright";
    CHECK(utility::extract_number<string>(s) == "92");
    CHECK(utility::extract_number<int>(s) == 92);

    s = "okay now it should be 311.51";
    CHECK(utility::extract_number<string>(s) == "311.51");
    CHECK(utility::extract_number<double>(s) == 311.51);

    s = "ad713e15.c";
    CHECK(utility::extract_number<string>(s) == "713");
    CHECK(utility::extract_number<int>(s) == 713);

    s = "814.98.txt";
    CHECK(utility::extract_number<string>(s) == "814.98");
    CHECK(utility::extract_number<double>(s) == 814.98);
}

TEST_CASE("utility_limits", "[utility]") {
    Limit lim1(0, 1);
    Limit lim2(-5, 15.5);

    SECTION("span") {
        CHECK(lim1.span() == 1);
        CHECK(lim2.span() == 20.5);
    }

    SECTION("center") {
        CHECK(lim1.center() == 0.5);
        CHECK(lim2.center() == 10.25);
    }

    SECTION("merge") {
        SECTION("simple") {
            lim1.merge(lim2);
            CHECK(lim1.min == -5);
            CHECK(lim1.max == 15.5);
        }

        SECTION("overlap") {
            Limit lim2(0.5, 1.5);
            lim1.merge(lim2);
            CHECK(lim1.min == 0);
            CHECK(lim1.max == 1.5);
        }
    }

    SECTION("expand") {
        lim1.expand(0.1);
        CHECK(lim1.min == -0.1);
        CHECK(lim1.max == 1.1);

        lim1 = Limit(0, 5);
        lim1.expand(0.1);
        CHECK(lim1.min == -0.5);
        CHECK(lim1.max == 5.5);
    }

    SECTION("add & subtract") {
        lim1 += 2;
        CHECK(lim1.min == 2);
        CHECK(lim1.max == 3);

        lim1 -= 0.5;
        CHECK(lim1.min == 1.5);
        CHECK(lim1.max == 2.5);
    }
}