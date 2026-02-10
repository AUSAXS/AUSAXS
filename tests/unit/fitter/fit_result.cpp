#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <fitter/FitResult.h>
#include <mini/detail/Result.h>
#include <mini/detail/FittedParameter.h>

using namespace ausaxs;
using namespace ausaxs::fitter;

TEST_CASE("FitResult::constructor") {
    SECTION("from mini::Result") {
        mini::Result res;
        res.fval = 1.5;
        res.status = 0;
        res.fevals = 100;
        
        FitResult result(res, 5);
        REQUIRE(result.fval == 1.5);
        REQUIRE(result.dof == 5);
        REQUIRE(result.status == 0);
        REQUIRE(result.fevals == 100);
    }
}

TEST_CASE("FitResult::add_fit") {
    FitResult result1;
    result1.dof = 10;
    result1.parameters = {{"a", 1.0, 0.1}, {"b", 2.0, 0.2}};
    
    FitResult result2;
    result2.parameters = {{"c", 3.0, 0.3}};
    
    SECTION("add to back") {
        result1.add_fit(&result2, false);
        REQUIRE(result1.parameters.size() == 3);
        REQUIRE(result1.parameters[0].name == "a");
        REQUIRE(result1.parameters[1].name == "b");
        REQUIRE(result1.parameters[2].name == "c");
        REQUIRE(result1.dof == 9);
    }

    SECTION("add to front") {
        result1.add_fit(&result2, true);
        REQUIRE(result1.parameters.size() == 3);
        REQUIRE(result1.parameters[0].name == "c");
        REQUIRE(result1.parameters[1].name == "a");
        REQUIRE(result1.parameters[2].name == "b");
        REQUIRE(result1.dof == 9);
    }
}

TEST_CASE("FitResult::set_data_curves") {
    FitResult result;
    
    SECTION("with vectors") {
        std::vector<double> q = {0.1, 0.2, 0.3};
        std::vector<double> data = {1.0, 2.0, 3.0};
        std::vector<double> data_err = {0.1, 0.1, 0.1};
        std::vector<double> model = {1.1, 2.1, 3.1};
        std::vector<double> residuals = {0.1, 0.1, 0.1};
        
        result.set_data_curves(
            std::vector<double>(q),
            std::vector<double>(data),
            std::vector<double>(data_err),
            std::vector<double>(model),
            std::vector<double>(residuals)
        );
        
        REQUIRE(result.curves.size_cols() == 5);
        REQUIRE(result.curves.is_named());
        auto names = result.curves.get_col_names();
        REQUIRE(names.size() == 5);
        REQUIRE(names[0] == "q");
        REQUIRE(names[1] == "I");
        REQUIRE(names[2] == "I_err");
        REQUIRE(names[3] == "I_fit");
        REQUIRE(names[4] == "residuals");
    }
}

TEST_CASE("FitResult::to_string") {
    FitResult result;
    result.status = 0;
    result.fevals = 100;
    result.fval = 10.5;
    result.dof = 5;
    result.parameters = {{"a", 1.23, 0.05}, {"b", 4.56, 0.10}};
    
    std::string str = result.to_string();
    
    REQUIRE(str.find("FIT REPORT") != std::string::npos);
    REQUIRE(str.find("Converged: yes") != std::string::npos);
    REQUIRE(str.find("Fevals:") != std::string::npos);
    REQUIRE(str.find("chi2:") != std::string::npos);
    REQUIRE(str.find("dof:") != std::string::npos);
    REQUIRE(str.find("a") != std::string::npos);
    REQUIRE(str.find("b") != std::string::npos);
}
