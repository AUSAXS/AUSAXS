#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <minimizer/ROOT.h>
#include <minimizer/Scan.h>
#include <minimizer/Golden.h>
#include <plots/all.h>

struct UnitaryTestFunction {
    UnitaryTestFunction(std::function<double(double*)> function, Limit bounds, double min) : function(function), bounds(bounds), min(min) {}

    std::function<double(double*)> function;
    Limit bounds;
    double min;
};


UnitaryTestFunction problem04([] (double* pars) {double x = pars[0]; return -(16*x*x - 24*x + 5)*std::exp(-x);}, Limit(1, 6), 2.868034);
UnitaryTestFunction problem13([] (double* pars) {double x = pars[0]; return -std::pow(x, 0.66) - std::pow(1 - x*x, 0.33);}, Limit(0, 1), 1./std::sqrt(2));
UnitaryTestFunction problem18([] (double* pars) {double x = pars[0]; return x <= 3 ? (x-2)*(x-2) : 2*std::log(x - 2) + 1;}, Limit(0, 6), 2);

TEST_CASE("golden_landscape", "[minimizer],[manual]") {
    mini::Golden mini(problem04.function, "a", problem04.bounds);
    auto res = mini.minimize();

    Dataset evaluations = mini.get_evaluated_points();
    Dataset landscape = mini.landscape();
    landscape.add_plot_options("lines", {{"color", kBlack}});
    evaluations.add_plot_options("markers", {{"color", kOrange+2}});

    plots::PlotDataset plot(landscape);
    plot.plot(evaluations);
    plot.save("figures/test/minimizer/golden_test.pdf");
}

TEST_CASE("golden_minimizer", "[minimizer]") {
    auto GoldenTest = [] (const UnitaryTestFunction& test) {
        mini::Golden mini(test.function, "a", test.bounds);
        auto res = mini.minimize();
        CHECK_THAT(res.parameters.at("a"), Catch::Matchers::WithinRel(test.min, mini.tol));
    };

    GoldenTest(problem04);
    GoldenTest(problem13);
    GoldenTest(problem18);
}

TEST_CASE("scan_minimizer", "[minimizer]") {
}

TEST_CASE("root_minimizer", "[minimizer]") {
}