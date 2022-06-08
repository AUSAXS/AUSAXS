#include <catch2/catch_test_macros.hpp>
#include <fitter/FitReporter.h>
#include <plots/PlotOptions.h>
#include <plots/PlotDataset.h>
#include <utility/Utility.h>

using std::vector, std::string;

TEST_CASE("plots", "[utility],[manual]") {
    Dataset data;
    data.x = {-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5};
    data.y = {-5, -3, 4, 8, 12, 7, 3, 1, -3, -5, -9};
    data.yerr = vector<double>(data.size(), 0.5);

    SECTION("markers") {
        plots::PlotDataset::quick_plot(data, "figures/test/utility/plots/default.png");

        data.add_plot_options("markers");
        plots::PlotDataset::quick_plot(data, "figures/test/utility/plots/markers.png");

        data.add_plot_options("errors");
        plots::PlotDataset::quick_plot(data, "figures/test/utility/plots/errors.png");
    }

    SECTION("limits") {
        data.add_plot_options("markers", {{"ylimits", Limit(3, 5)}});
        plots::PlotDataset::quick_plot(data, "figures/test/utility/plots/limits.png");
    }

    SECTION("log") {
        data.ylimits(1, inf);
        data.add_plot_options({{"logy", true}});
        plots::PlotDataset::quick_plot(data, "figures/test/utility/plots/log.png");
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

TEST_CASE("plotoptions", "[utility]") {
    plots::PlotOptions options;

    SECTION("set") {
        options.set({{"xlabel", "label for x"}, {"ylabel", string("label for y")}, {"drawline", true}, {"draw_errors", true}, {"color", kBlue}, {"alpha", 0.3}});
        CHECK(options.xlabel == "label for x");
        CHECK(options.ylabel == "label for y");
        CHECK(options.draw_line == true);
        CHECK(options.draw_errors == true);
        CHECK(options.color == kBlue);
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