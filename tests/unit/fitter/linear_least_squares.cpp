#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <fitter/detail/LinearLeastSquares.h>

using namespace ausaxs;
using namespace ausaxs::fitter::detail;

TEST_CASE("LinearLeastSquares::constructor") {
    SECTION("unity errors") {
        std::vector<double> data = {1.0, 2.0, 3.0};
        std::vector<double> model = {2.0, 4.0, 6.0};
        LinearLeastSquares fitter(data, model);
        
        REQUIRE(fitter.size() == 3);
        REQUIRE(fitter.dof() == 1);
    }

    SECTION("custom errors") {
        std::vector<double> data = {1.0, 2.0, 3.0};
        std::vector<double> model = {2.0, 4.0, 6.0};
        std::vector<double> errors = {0.1, 0.1, 0.1};
        LinearLeastSquares fitter(data, model, errors);
        
        REQUIRE(fitter.size() == 3);
        REQUIRE(fitter.dof() == 1);
    }
}

TEST_CASE("LinearLeastSquares::size") {
    std::vector<double> data = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> model = {2.0, 4.0, 6.0, 8.0, 10.0};
    LinearLeastSquares fitter(data, model);
    
    REQUIRE(fitter.size() == 5);
}

TEST_CASE("LinearLeastSquares::dof") {
    SECTION("3 points") {
        std::vector<double> data = {1.0, 2.0, 3.0};
        std::vector<double> model = {2.0, 4.0, 6.0};
        LinearLeastSquares fitter(data, model);
        
        REQUIRE(fitter.dof() == 1);
    }

    SECTION("5 points") {
        std::vector<double> data = {1.0, 2.0, 3.0, 4.0, 5.0};
        std::vector<double> model = {2.0, 4.0, 6.0, 8.0, 10.0};
        LinearLeastSquares fitter(data, model);
        
        REQUIRE(fitter.dof() == 3);
    }

    SECTION("10 points") {
        std::vector<double> data(10, 1.0);
        std::vector<double> model(10, 2.0);
        LinearLeastSquares fitter(data, model);
        
        REQUIRE(fitter.dof() == 8);
    }
}

TEST_CASE("LinearLeastSquares::fit") {
    SECTION("perfect linear relationship y = 2x") {
        std::vector<double> data = {1.0, 2.0, 3.0, 4.0};
        std::vector<double> model = {2.0, 4.0, 6.0, 8.0};
        LinearLeastSquares fitter(data, model);
        
        auto result = fitter.fit();
        REQUIRE(result != nullptr);
        REQUIRE(result->dof == 2);
        REQUIRE(result->parameters.size() == 2);
        REQUIRE(result->parameters[0].name == "a");
        REQUIRE(result->parameters[1].name == "b");
        REQUIRE_THAT(result->parameters[0].value, Catch::Matchers::WithinAbs(2.0, 1e-10));
        REQUIRE_THAT(result->parameters[1].value, Catch::Matchers::WithinAbs(0.0, 1e-10));
    }

    SECTION("linear relationship y = 2x + 1") {
        std::vector<double> data = {1.0, 2.0, 3.0, 4.0};
        std::vector<double> model = {3.0, 5.0, 7.0, 9.0};
        LinearLeastSquares fitter(data, model);
        
        auto result = fitter.fit();
        REQUIRE(result->parameters.size() == 2);
        REQUIRE_THAT(result->parameters[0].value, Catch::Matchers::WithinAbs(2.0, 1e-10));
        REQUIRE_THAT(result->parameters[1].value, Catch::Matchers::WithinAbs(1.0, 1e-10));
    }

    SECTION("with errors") {
        std::vector<double> data = {1.0, 2.0, 3.0, 4.0};
        std::vector<double> model = {2.0, 4.0, 6.0, 8.0};
        std::vector<double> errors = {0.1, 0.1, 0.1, 0.1};
        LinearLeastSquares fitter(data, model, errors);
        
        auto result = fitter.fit();
        REQUIRE(result->parameters.size() == 2);
        REQUIRE_THAT(result->parameters[0].value, Catch::Matchers::WithinAbs(2.0, 1e-10));
        REQUIRE_THAT(result->parameters[1].value, Catch::Matchers::WithinAbs(0.0, 1e-10));
    }
}

TEST_CASE("LinearLeastSquares::get_residuals") {
    std::vector<double> data = {1.0, 2.0, 3.0};
    std::vector<double> model = {2.0, 4.0, 6.0};
    LinearLeastSquares fitter(data, model);
    
    SECTION("perfect fit") {
        auto residuals = fitter.get_residuals({2.0, 0.0});
        REQUIRE(residuals.size() == 3);
        for (auto r : residuals) {
            REQUIRE_THAT(r, Catch::Matchers::WithinAbs(0.0, 1e-10));
        }
    }

    SECTION("imperfect fit") {
        auto residuals = fitter.get_residuals({1.5, 0.5});
        REQUIRE(residuals.size() == 3);
        REQUIRE_THAT(residuals[0], Catch::Matchers::WithinAbs(2.0 - 2.0, 1e-10));
        REQUIRE_THAT(residuals[1], Catch::Matchers::WithinAbs(4.0 - 3.5, 1e-10));
        REQUIRE_THAT(residuals[2], Catch::Matchers::WithinAbs(6.0 - 5.0, 1e-10));
    }
}

TEST_CASE("LinearLeastSquares::get_model_curve") {
    std::vector<double> data = {1.0, 2.0, 3.0};
    std::vector<double> model = {2.0, 4.0, 6.0};
    LinearLeastSquares fitter(data, model);
    
    SECTION("with parameters") {
        auto curve = fitter.get_model_curve({2.0, 1.0});
        REQUIRE(curve.size() == 3);
        REQUIRE_THAT(curve[0], Catch::Matchers::WithinAbs(3.0, 1e-10));
        REQUIRE_THAT(curve[1], Catch::Matchers::WithinAbs(5.0, 1e-10));
        REQUIRE_THAT(curve[2], Catch::Matchers::WithinAbs(7.0, 1e-10));
    }

    SECTION("without parameters") {
        auto curve = fitter.get_model_curve();
        REQUIRE(curve.size() == 3);
    }
}



TEST_CASE("LinearLeastSquares::chi2") {
    std::vector<double> data = {1.0, 2.0, 3.0};
    std::vector<double> model = {2.0, 4.0, 6.0};
    LinearLeastSquares fitter(data, model);
    
    SECTION("perfect fit") {
        double chi2 = fitter.chi2({2.0, 0.0});
        REQUIRE_THAT(chi2, Catch::Matchers::WithinAbs(0.0, 1e-10));
    }

    SECTION("imperfect fit") {
        double chi2 = fitter.chi2({1.0, 0.0});
        REQUIRE(chi2 > 0.0);
    }
}
