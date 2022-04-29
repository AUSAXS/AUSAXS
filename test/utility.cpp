#include <catch2/catch_test_macros.hpp>
#include <fitter/FitReporter.h>
#include <plots/PlotOptions.h>

using std::vector, std::string;

TEST_CASE("fitreporter", "[utility],[manual]") {
    Fit fit;
    fit.converged = true;
    fit.calls = 100;
    fit.chi2 = 1000;
    fit.dof = 3;
    fit.params = {{"a", 2}, {"b", 3}, {"c", 4}};
    fit.errors = {{"a", 0.5}, {"b", 0.4}, {"c", 0.3}};

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
        fit.converged = false;
        fit2.calls = 20;
        fit2.chi2 = 200;
        fit2.dof = 6;
        fit2.params = {{"a", 1}, {"b", 0.5}, {"c", 8}};
        fit2.errors = {{"a", 0.1}, {"b", 0.3}, {"c", 0.6}};

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