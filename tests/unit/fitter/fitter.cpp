#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <fitter/Fitter.h>
#include <mini/MiniFwd.h>

using namespace ausaxs;
using namespace ausaxs::fitter;

class TestFitter : public Fitter {
public:
    TestFitter(std::vector<double> data, std::vector<double> model) 
        : data(data), model(model) {}
    
    std::unique_ptr<FitResult> fit() override {
        return std::make_unique<FitResult>();
    }
    
    std::vector<double> fit_params_only() override {
        return {1.0, 0.0};
    }
    
    unsigned int dof() const override {
        return data.size() - 2;
    }
    
    unsigned int size() const override {
        return data.size();
    }
    
    std::vector<double> get_residuals(const std::vector<double>& params) override {
        std::vector<double> residuals(data.size());
        for (size_t i = 0; i < data.size(); ++i) {
            residuals[i] = model[i] - (params[0] * data[i] + params[1]);
        }
        return residuals;
    }
    
private:
    std::vector<double> data;
    std::vector<double> model;
};

TEST_CASE("Fitter::dof") {
    std::vector<double> data = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> model = {2.0, 4.0, 6.0, 8.0, 10.0};
    TestFitter fitter(data, model);
    
    REQUIRE(fitter.dof() == 3);
}

TEST_CASE("Fitter::degrees_of_freedom") {
    std::vector<double> data = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> model = {2.0, 4.0, 6.0, 8.0, 10.0};
    TestFitter fitter(data, model);
    
    REQUIRE(fitter.degrees_of_freedom() == fitter.dof());
    REQUIRE(fitter.degrees_of_freedom() == 3);
}

TEST_CASE("Fitter::size") {
    SECTION("5 points") {
        std::vector<double> data = {1.0, 2.0, 3.0, 4.0, 5.0};
        std::vector<double> model = {2.0, 4.0, 6.0, 8.0, 10.0};
        TestFitter fitter(data, model);
        
        REQUIRE(fitter.size() == 5);
    }
    
    SECTION("3 points") {
        std::vector<double> data = {1.0, 2.0, 3.0};
        std::vector<double> model = {2.0, 4.0, 6.0};
        TestFitter fitter(data, model);
        
        REQUIRE(fitter.size() == 3);
    }
}

TEST_CASE("Fitter::chi2") {
    std::vector<double> data = {1.0, 2.0, 3.0};
    std::vector<double> model = {2.0, 4.0, 6.0};
    TestFitter fitter(data, model);
    
    SECTION("perfect fit") {
        double chi2 = fitter.chi2({2.0, 0.0});
        REQUIRE_THAT(chi2, Catch::Matchers::WithinAbs(0.0, 1e-10));
    }
    
    SECTION("imperfect fit") {
        double chi2 = fitter.chi2({1.5, 0.5});
        REQUIRE(chi2 > 0.0);
        
        double expected = 0.0;
        for (size_t i = 0; i < data.size(); ++i) {
            double residual = model[i] - (1.5 * data[i] + 0.5);
            expected += residual * residual;
        }
        REQUIRE_THAT(chi2, Catch::Matchers::WithinAbs(expected, 1e-10));
    }
}

TEST_CASE("Fitter::fit_chi2_only") {
    std::vector<double> data = {1.0, 2.0, 3.0};
    std::vector<double> model = {1.0, 2.0, 3.0};
    TestFitter fitter(data, model);
    
    double chi2 = fitter.fit_chi2_only();
    REQUIRE(chi2 >= 0.0);
}

TEST_CASE("Fitter::set_algorithm") {
    std::vector<double> data = {1.0, 2.0, 3.0};
    std::vector<double> model = {2.0, 4.0, 6.0};
    TestFitter fitter(data, model);
    
    fitter.set_algorithm(mini::algorithm::BFGS);
}
