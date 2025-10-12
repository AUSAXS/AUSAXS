#include <catch2/catch_test_macros.hpp>
#include <fitter/FitReporter.h>
#include <plots/PlotOptions.h>
#include <plots/PlotDataset.h>
#include <utility/Utility.h>
#include <utility/StringUtils.h>
#include <settings/GeneralSettings.h>
#include <plots/All.h>

using namespace ausaxs;

TEST_CASE("fitreporter", "[manual]") {
    fitter::FitResult fit;
    fit.status = true;
    fit.fevals = 100;
    fit.fval = 1000;
    fit.dof = 3;
    fit.parameters = {{"a", 2, 0.5}, {"b", 3, 0.4}, {"c", 4, 0.3}};

    SECTION("Single") {
        SECTION("Terminal printing") {
            fitter::FitReporter::report(&fit);
        }
        SECTION("File printing") {
            fitter::FitReporter::save(&fit, "temp/fitreport1.txt");
        }
    }

    SECTION("Multi") {
        fitter::FitResult fit2;
        fit.status = false;
        fit2.fevals = 20;
        fit2.fval = 200;
        fit2.dof = 6;
        fit2.parameters = {{"a", 1, 0.1}, {"b", 0.5, 0.3}, {"c", 8, 0.6}};

        std::vector<fitter::FitResult> fits = {fit, fit2};
        std::vector<std::string> titles = {"First fit", "Second fit"};

        SECTION("Terminal printing") {
            fitter::FitReporter::report(fits, titles);   
        }
        SECTION("File printing") {
            fitter::FitReporter::save(fits, "temp/fitreport2.txt", titles);
        }
    }
}

TEST_CASE("plotoptions") {
    plots::PlotOptions options;

    SECTION("set") {
        options.set({{"xlabel", "label for x"}, {"ylabel", std::string("label for y")}, {"line", true}, {"errors", true}, {"color", style::color::blue}, {"alpha", 0.3}});
        CHECK(options.xlabel == "label for x");
        CHECK(options.ylabel == "label for y");
        CHECK(options.draw_line == true);
        CHECK(options.draw_errors == true);
        CHECK(options.color == style::color::blue);
        CHECK(options.alpha == 0.3);
    }
}

TEST_CASE("limits") {  
    Limit lim1(0, 1);
    Limit lim2(-5, 15.5);

    SECTION("span") {
        CHECK(lim1.span() == 1);
        CHECK(lim2.span() == 20.5);
    }

    SECTION("center") {
        CHECK(lim1.center() == 0.5);
        CHECK(lim2.center() == 5.25);
    }

    SECTION("merge") {
        SECTION("simple") {
            lim1.merge(lim2);
            CHECK(lim1.min == -5);
            CHECK(lim1.max == 15.5);
        }

        SECTION("overlap") {
            Limit lim3(0.5, 1.5);
            lim1.merge(lim3);
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