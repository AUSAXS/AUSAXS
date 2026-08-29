#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <fitter/FitReporter.h>
#include <fitter/FitResult.h>
#include <mini/detail/FittedParameter.h>

#include <support/temp_file.h>

#include <fstream>

using namespace ausaxs;
using namespace ausaxs::fitter;

TEST_CASE("FitReporter::save") {
    FitResult result;
    result.status = 0;
    result.fevals = 50;
    result.fval = 5.5;
    result.dof = 3;
    result.parameters = {{"param1", 1.0, 0.1}};
    
    SECTION("with header") {
        test::TempFile path(".txt");
        FitReporter::save(&result, path, "Test Header");

        REQUIRE(path.exists());

        std::ifstream file(path.str());
        std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());

        REQUIRE_THAT(content, Catch::Matchers::ContainsSubstring("Test Header"));
        REQUIRE_THAT(content, Catch::Matchers::ContainsSubstring("FIT REPORT"));
    }

    SECTION("without header") {
        test::TempFile path(".txt");
        FitReporter::save(&result, path);

        REQUIRE(path.exists());

        std::ifstream file(path.str());
        std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());

        REQUIRE_THAT(content, Catch::Matchers::ContainsSubstring("FIT REPORT"));
    }
}

TEST_CASE("FitReporter::save_multiple") {
    FitResult result1;
    result1.status = 0;
    result1.fevals = 50;
    result1.fval = 5.5;
    result1.dof = 3;
    result1.parameters = {{"param1", 1.0, 0.1}};
    
    FitResult result2;
    result2.status = 0;
    result2.fevals = 75;
    result2.fval = 3.2;
    result2.dof = 4;
    result2.parameters = {{"param2", 2.0, 0.2}};
    
    SECTION("with titles") {
        test::TempFile path(".txt");
        std::vector<FitResult> fits = {result1, result2};
        std::vector<std::string> titles = {"Fit 1", "Fit 2"};

        FitReporter::save(fits, path, titles);

        REQUIRE(path.exists());

        std::ifstream file(path.str());
        std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());

        REQUIRE_THAT(content, Catch::Matchers::ContainsSubstring("Fit 1"));
        REQUIRE_THAT(content, Catch::Matchers::ContainsSubstring("Fit 2"));
    }

    SECTION("without titles") {
        test::TempFile path(".txt");
        std::vector<FitResult> fits = {result1, result2};

        FitReporter::save(fits, path);

        REQUIRE(path.exists());

        std::ifstream file(path.str());
        std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());

        REQUIRE_THAT(content, Catch::Matchers::ContainsSubstring("FIT REPORT"));
    }
}
