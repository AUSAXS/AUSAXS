#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <minimizer/ROOT.h>
#include <minimizer/Scan.h>
#include <minimizer/Golden.h>
#include <plots/all.h>

struct TestFunction {
    TestFunction(std::function<double(double*)> function, std::vector<Limit> bounds, std::vector<double> min) : function(function), bounds(bounds), min(min) {}
    TestFunction(std::function<double(double*)> function, Limit bounds, double min) : TestFunction(function, vector{bounds}, vector{min}) {}

    std::function<double(double*)> function;
    std::vector<Limit> bounds;
    std::vector<double> min;
};

// 1D functions
TestFunction problem04([] (double* pars) {double x = pars[0]; return -(16*x*x - 24*x + 5)*std::exp(-x);}, Limit(1, 6), 2.868034);
TestFunction problem13([] (double* pars) {double x = pars[0]; return -std::pow(x, 0.66) - std::pow(1 - x*x, 0.33);}, Limit(0, 1), 1./std::sqrt(2));
TestFunction problem18([] (double* pars) {double x = pars[0]; return x <= 3 ? (x-2)*(x-2) : 2*std::log(x - 2) + 1;}, Limit(0, 6), 2);

// 2D functions (nice)
TestFunction Decanomial([] (double* pars) {
    double x1 = pars[0], x2 = pars[1]; 
    return 0.001*std::pow(std::abs(std::pow(x2, 4) + 12*std::pow(x2, 3) + 54*std::pow(x2, 2) + 108*x2 + 81) 
                        + std::abs(std::pow(x1, 10) - 20*std::pow(x1, 9) + 180*std::pow(x1, 8) - 960*std::pow(x1, 7) + 3360*std::pow(x1, 6) - 
                              8064*std::pow(x1, 5) + 13340*std::pow(x1, 4) - 15360*std::pow(x1, 3) + 11520*std::pow(x1, 2) - 5120*x1 + 2624), 2);}, 
    {Limit(0, 2.5), Limit(-4, -2)}, 
    {2, -3}
);
TestFunction Hosaki([] (double* pars) {
    double x1 = pars[0], x2 = pars[1]; 
    return (1 - 8*x1 + 7*x1*x1 - 7*std::pow(x1, 3)/3 + std::pow(x1, 4)/4)*x2*x2*std::exp(-x1);}, 
    {Limit(0, 5), Limit(0, 5)}, 
    {4, 2}
);
TestFunction RosenbrockModified([] (double* pars) {
    double x1 = pars[0], x2 = pars[1]; 
    return 74 + 100*std::pow(x2 - x1*x1, 2) + std::pow(1-x1, 2) - 400*std::exp(-(std::pow(x1+1, 2)+std::pow(x2+1, 2))*10);}, 
    {Limit(-1, 0.4), Limit(-1, 1)}, 
    {-0.9, -0.95}
);

TEST_CASE("1d_landscape", "[minimizer],[manual]") {
    mini::Golden mini(problem04.function, {"x1", problem04.bounds[0]});
    auto res = mini.minimize();

    Dataset evaluations = mini.get_evaluated_points();
    Dataset landscape = mini.landscape();
    landscape.add_plot_options("lines", {{"color", kBlack}});
    evaluations.add_plot_options("markers", {{"color", kOrange+2}});

    plots::PlotDataset plot(landscape);
    plot.plot(evaluations);
    plot.save("figures/test/minimizer/golden_test.pdf");
}

// TEST_CASE("2d_landscape", "[minimizer],[manual]") {
//     mini::ROOT mini(Decanomial.function, {"x1", "x2"}, Decanomial.bounds);
//     auto res = mini.minimize();
// }

TEST_CASE("golden_minimizer", "[minimizer]") {
    auto GoldenTest = [] (const TestFunction& test) {
        mini::Golden mini(test.function, {"a", test.bounds[0]});
        auto res = mini.minimize();
        CHECK_THAT(res.get_parameter("a"), Catch::Matchers::WithinRel(test.min[0], mini.tol));
    };

    GoldenTest(problem04);
    GoldenTest(problem13);
    GoldenTest(problem18);
}

TEST_CASE("scan_minimizer", "[minimizer]") {
}

TEST_CASE("root_minimizer", "[minimizer]") {
    auto ROOTTest = [] (const TestFunction& test, std::pair<double, double> guess) {
        mini::ROOT mini(test.function, {{"x1", guess.first}, {"x2", guess.second}}, test.bounds);
        auto res = mini.minimize();
        CHECK_THAT(res.get_parameter("x1"), Catch::Matchers::WithinRel(test.min[0], mini.tol));
        CHECK_THAT(res.get_parameter("x2"), Catch::Matchers::WithinRel(test.min[1], mini.tol));
    };

    GoldenTest(Decanomial);
    GoldenTest(Hosaki);
    GoldenTest(RosenbrockModified);
}