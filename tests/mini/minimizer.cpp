#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <mini/detail/Parameter.h>

#include <mini/All.h>
#include <plots/All.h>

using std::vector;
using namespace ausaxs;

struct TestFunction {
    TestFunction(std::function<double(std::vector<double>)> function, const std::vector<Limit>& bounds, const std::vector<double>& min) : function(function), bounds(bounds), min(min) {}
    TestFunction(std::function<double(std::vector<double>)> function, const Limit& bounds, double min) : TestFunction(function, vector{bounds}, vector{min}) {}

    std::vector<double> get_center() const {
        std::vector<double> v;
        std::for_each(bounds.begin(), bounds.end(), [&v] (const Limit& lim) {v.push_back(lim.center());});
        return v;
    }

    std::function<double(std::vector<double>)> function;
    std::vector<Limit> bounds;
    std::vector<double> min;
};

// 1D functions
TestFunction problem04([] (std::vector<double> pars) {double x = pars[0]; return -(16*x*x - 24*x + 5)*std::exp(-x);}, Limit(1, 6), 2.868034);
TestFunction problem13([] (std::vector<double> pars) {double x = pars[0]; return -std::pow(x, 0.66) - std::pow(1 - x*x, 0.33);}, Limit(0, 1), 1./std::sqrt(2));
TestFunction problem18([] (std::vector<double> pars) {double x = pars[0]; return x <= 3 ? (x-2)*(x-2) : 2*std::log(x - 2) + 1;}, Limit(0, 6), 2);

// 2D functions (nice)
TestFunction Decanomial([] (std::vector<double> pars) {
    double x1 = pars[0], x2 = pars[1]; 
    return 0.001*std::pow(std::abs(std::pow(x2, 4) + 12*std::pow(x2, 3) + 54*std::pow(x2, 2) + 108*x2 + 81) 
                        + std::abs(std::pow(x1, 10) - 20*std::pow(x1, 9) + 180*std::pow(x1, 8) - 960*std::pow(x1, 7) + 3360*std::pow(x1, 6) - 
                              8064*std::pow(x1, 5) + 13340*std::pow(x1, 4) - 15360*std::pow(x1, 3) + 11520*std::pow(x1, 2) - 5120*x1 + 2624), 2);}, 
    {Limit(0, 2.5), Limit(-4, -2)}, 
    {2, -3}
);
TestFunction Hosaki([] (std::vector<double> pars) {
    double x1 = pars[0], x2 = pars[1]; 
    return (1 - 8*x1 + 7*x1*x1 - 7*std::pow(x1, 3)/3 + std::pow(x1, 4)/4)*x2*x2*std::exp(-x2);}, 
    {Limit(0, 5), Limit(0, 5)}, 
    {4, 2}
);
TestFunction Rosenbrock([] (std::vector<double> pars) {
    double x1 = pars[0], x2 = pars[1]; 
    return std::pow(1-x1, 2) + 100*std::pow(x2-x1*x1, 2);}, 
    {Limit(-2, 2), Limit(-2, 2)}, 
    {1, 1}
);

TEST_CASE("1d_landscape", "[manual]") {
    mini::Golden mini(problem04.function, {"x1", problem04.bounds[0]});
    auto res = mini.minimize();

    SimpleDataset evaluations = mini.get_evaluated_points().as_dataset();
    SimpleDataset landscape = mini.landscape().as_dataset();

    plots::PlotDataset plot(landscape, plots::PlotOptions(style::draw::line, {{"color", style::color::black}}));
    plot.plot(evaluations, plots::PlotOptions(style::draw::points, {{"color", style::color::orange}}));
    plot.save("figures/tests/minimizer/golden_test.pdf");
}

TEST_CASE("Minimizer: golden") {
    auto GoldenTest = [] (const TestFunction& test) {
        mini::Golden mini(test.function, {"a", test.bounds[0]});
        auto res = mini.minimize();
        CHECK_THAT(res.get_parameter("a").value, Catch::Matchers::WithinAbs(test.min[0], mini.tol));
    };

    SECTION("problem04") {GoldenTest(problem04);}
    SECTION("problem13") {GoldenTest(problem13);}
    SECTION("problem18") {GoldenTest(problem18);}
}

TEST_CASE("Minimizer: scan") {
    auto ScanTest1D = [] (const TestFunction& test) {
        mini::Scan mini(test.function, {"a", test.bounds[0]});
        mini.set_max_evals(1000);
        auto res = mini.minimize();
        CHECK_THAT(res.get_parameter("a").value, Catch::Matchers::WithinAbs(test.min[0], mini.tol));
    };

    auto ScanTest1DRough = [] (const TestFunction& test) {
        mini::Scan mini(test.function, {"a", test.bounds[0]});
        mini.set_max_evals(10);
        auto res = mini.minimize();
        CHECK_THAT(res.get_parameter("a").value, Catch::Matchers::WithinAbs(test.min[0], mini.tol));
    };

    // test with a fine grid
    SECTION("problem04") {ScanTest1D(problem04);}
    SECTION("problem13") {ScanTest1D(problem13);}
    SECTION("problem18") {ScanTest1D(problem18);}

    // test with a rough grid & let the local minimizer find the actual minima
    SECTION("problem04 rough") {ScanTest1DRough(problem04);}
    SECTION("problem13 rough") {ScanTest1DRough(problem13);}
    SECTION("problem18 rough") {ScanTest1DRough(problem18);}
}

// TEST_CASE("Minimizer: minimum_explorer") {
//     auto ExplorerTest1D = [] (const TestFunction& test) {
//         mini::dlibMinimizer<mini::algorithm::BFGS> mini1(test.function, {{"a", test.bounds[0]}});
//         auto res = mini1.minimize();

//         mini::Parameter p = res.get_parameter("a");
//         mini::MinimumExplorer mini2(test.function, p, 100);
//         res = mini2.minimize();
//         mini::Golden mini3(test.function, {"a", test.bounds[0]});
//         SimpleDataset line = mini3.landscape(1000).as_dataset();
//         CHECK_THAT(res.get_parameter("a").value, Catch::Matchers::WithinAbs(test.min[0], mini1.tol));
//     };

//     // test with a fine grid
//     SECTION("problem04") {ExplorerTest1D(problem04);}
//     SECTION("problem13") {ExplorerTest1D(problem13);}
//     SECTION("problem18") {ExplorerTest1D(problem18);}
// }

#ifdef DLIB_AVAILABLE
TEST_CASE("Minimizer: dlib") {
    auto dlibTest1D = [] (const TestFunction& test, mini::algorithm type) {
        if (type == mini::algorithm::BFGS) {
            auto mini = mini::dlibMinimizer<mini::algorithm::BFGS>(test.function, {mini::Parameter{"a", test.get_center()[0], test.bounds[0]}});
            auto res = mini.minimize();
            CHECK_THAT(res.get_parameter("a").value, Catch::Matchers::WithinAbs(test.min[0], mini.tol));
        } else if (type == mini::algorithm::DLIB_GLOBAL) {
            auto mini = mini::dlibMinimizer<mini::algorithm::DLIB_GLOBAL>(test.function, {mini::Parameter{"a", test.get_center()[0], test.bounds[0]}});
            auto res = mini.minimize();
            CHECK_THAT(res.get_parameter("a").value, Catch::Matchers::WithinAbs(test.min[0], mini.tol));
        }
    };

    auto dlibTest2D = [] (const TestFunction& test, mini::algorithm type) {
        if (type == mini::algorithm::BFGS) {
            auto mini = mini::dlibMinimizer<mini::algorithm::BFGS>(test.function, {mini::Parameter{"a", test.bounds[0].center(), test.bounds[0]}, mini::Parameter{"b", test.bounds[1].center(), test.bounds[1]}});
            auto res = mini.minimize();
            CHECK_THAT(res.get_parameter("a").value, Catch::Matchers::WithinAbs(test.min[0], mini.tol));
            CHECK_THAT(res.get_parameter("b").value, Catch::Matchers::WithinAbs(test.min[1], mini.tol));
        } else if (type == mini::algorithm::DLIB_GLOBAL) {
            auto mini = mini::dlibMinimizer<mini::algorithm::BFGS>(test.function, {mini::Parameter{"a", test.bounds[0].center(), test.bounds[0]}, mini::Parameter{"b", test.bounds[1].center(), test.bounds[1]}});
            auto res = mini.minimize();
            CHECK_THAT(res.get_parameter("a").value, Catch::Matchers::WithinAbs(test.min[0], mini.tol));
            CHECK_THAT(res.get_parameter("b").value, Catch::Matchers::WithinAbs(test.min[1], mini.tol));
        }
    };

    SECTION("bfgs") {
        SECTION("problem04") {dlibTest1D(problem04, mini::algorithm::BFGS);}
        SECTION("problem13") {dlibTest1D(problem13, mini::algorithm::BFGS);}
        SECTION("problem18") {dlibTest1D(problem18, mini::algorithm::BFGS);}

        SECTION("Decanomial") {dlibTest2D(Decanomial, mini::algorithm::BFGS);}
        SECTION("Hosaki")     {dlibTest2D(Hosaki, mini::algorithm::BFGS);}
        SECTION("Rosenbrock") {dlibTest2D(Rosenbrock, mini::algorithm::BFGS);}
    }

    SECTION("dlib_global") {
        SECTION("problem04") {dlibTest1D(problem04, mini::algorithm::DLIB_GLOBAL);}
        SECTION("problem13") {dlibTest1D(problem13, mini::algorithm::DLIB_GLOBAL);}
        SECTION("problem18") {dlibTest1D(problem18, mini::algorithm::DLIB_GLOBAL);}

        SECTION("Decanomial") {dlibTest2D(Decanomial, mini::algorithm::DLIB_GLOBAL);}
        SECTION("Hosaki")     {dlibTest2D(Hosaki, mini::algorithm::DLIB_GLOBAL);}
        SECTION("Rosenbrock") {dlibTest2D(Rosenbrock, mini::algorithm::DLIB_GLOBAL);}
    }
}
#endif

TEST_CASE("Minimizer: create_minimizer") {
    #ifdef DLIB_AVAILABLE
        SECTION("dlib") {
            auto dlib = mini::create_minimizer(mini::algorithm::BFGS, problem04.function, {"a", problem04.bounds[0]});
            auto res = dlib->minimize();
            CHECK_THAT(res.get_parameter("a").value, Catch::Matchers::WithinAbs(problem04.min[0], dlib->tol));
        }
    #endif

    SECTION("golden") {
        auto golden = mini::create_minimizer(mini::algorithm::GOLDEN, problem04.function, {"a", problem04.bounds[0]});
        auto res = golden->minimize();
        CHECK_THAT(res.get_parameter("a").value, Catch::Matchers::WithinAbs(problem04.min[0], golden->tol));
    }
}