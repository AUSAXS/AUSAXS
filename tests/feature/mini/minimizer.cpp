#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <utility>

#include <mini/All.h>
#include <mini/detail/Parameter.h>
#include <numbers>
#include <plots/All.h>


using std::vector;
using namespace ausaxs;

// every minimized function is a sum of squared residuals, so the test problems are least-squares problems
struct TestFunction {
    TestFunction(mini::Minimizer::residual_function function, const std::vector<Limit>& bounds, const std::vector<double>& min) : function(std::move(function)), bounds(bounds), min(min) {}
    TestFunction(mini::Minimizer::residual_function function, const Limit& bounds, double min) : TestFunction(std::move(function), vector{bounds}, vector{min}) {}

    std::vector<double> get_center() const {
        std::vector<double> v;
        std::ranges::for_each(bounds, [&v] (const Limit& lim) {v.push_back(lim.center());});
        return v;
    }

    std::vector<mini::Parameter> get_parameters() const {
        std::vector<mini::Parameter> p;
        for (int i = 0; i < static_cast<int>(bounds.size()); ++i) {p.emplace_back("p" + std::to_string(i), bounds[i].center(), bounds[i]);}
        return p;
    }

    mini::Minimizer::residual_function function;
    std::vector<Limit> bounds;
    std::vector<double> min;
};

// noise-free samples of y = A exp(-k t), with A = 3 and k = 0.7
static std::vector<double> decay_data(double A, double k) {
    std::vector<double> y;
    for (int i = 0; i < 50; ++i) {y.push_back(A*std::exp(-k*0.1*i));}
    return y;
}
static std::vector<double> decay_residuals(double A, double k) {
    static const auto y = decay_data(3, 0.7);
    auto model = decay_data(A, k);
    for (int i = 0; i < static_cast<int>(y.size()); ++i) {model[i] -= y[i];}
    return model;
}

// 1D functions
static TestFunction decay1d([] (const std::vector<double>& p) {return decay_residuals(3, p[0]);}, Limit(0, 5), 0.7);
static TestFunction sqrt2([] (const std::vector<double>& p) {return std::vector{p[0]*p[0] - 2};}, Limit(0, 3), std::numbers::sqrt2);
static TestFunction euler([] (const std::vector<double>& p) {return std::vector{std::log(p[0]) - 1};}, Limit(0.5, 6), std::numbers::e);

// 2D functions
static TestFunction Rosenbrock([] (const std::vector<double>& p) {return std::vector{1 - p[0], 10*(p[1] - p[0]*p[0])};}, {Limit(-2, 2), Limit(-2, 2)}, {1, 1});
static TestFunction Beale([] (const std::vector<double>& p) {
    double x = p[0], y = p[1];
    return std::vector{1.5 - x + x*y, 2.25 - x + x*y*y, 2.625 - x + x*y*y*y};},
    {Limit(0, 4), Limit(0, 1)},
    {3, 0.5}
);
static TestFunction decay2d([] (const std::vector<double>& p) {return decay_residuals(p[0], p[1]);}, {Limit(0, 10), Limit(0, 5)}, {3, 0.7});

TEST_CASE("1d_landscape", "[manual]") {
    mini::Golden mini(decay1d.function, {"x1", decay1d.bounds[0]});
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

    SECTION("decay") {GoldenTest(decay1d);}
    SECTION("sqrt2") {GoldenTest(sqrt2);}
    SECTION("euler") {GoldenTest(euler);}
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
    SECTION("decay") {ScanTest1D(decay1d);}
    SECTION("sqrt2") {ScanTest1D(sqrt2);}
    SECTION("euler") {ScanTest1D(euler);}

    // test with a rough grid & let the local minimizer find the actual minima
    SECTION("decay rough") {ScanTest1DRough(decay1d);}
    SECTION("sqrt2 rough") {ScanTest1DRough(sqrt2);}
    SECTION("euler rough") {ScanTest1DRough(euler);}
}

#ifdef DLIB_AVAILABLE
TEST_CASE("Minimizer: dlib") {
    auto dlibTest = [] (const TestFunction& test, mini::algorithm type) {
        std::unique_ptr<mini::Minimizer> mini;
        if (type == mini::algorithm::BFGS) {
            mini = std::make_unique<mini::dlibMinimizer<mini::algorithm::BFGS>>(test.function, test.get_parameters());
        } else {
            mini = std::make_unique<mini::dlibMinimizer<mini::algorithm::DLIB_GLOBAL>>(test.function, test.get_parameters());
        }
        auto res = mini->minimize();
        for (int i = 0; i < static_cast<int>(test.min.size()); ++i) {
            CHECK_THAT(res.get_parameter(i).value, Catch::Matchers::WithinAbs(test.min[i], mini->tol));
        }
    };

    SECTION("bfgs") {
        SECTION("decay") {dlibTest(decay1d, mini::algorithm::BFGS);}
        SECTION("sqrt2") {dlibTest(sqrt2, mini::algorithm::BFGS);}
        SECTION("euler") {dlibTest(euler, mini::algorithm::BFGS);}

        SECTION("Rosenbrock") {dlibTest(Rosenbrock, mini::algorithm::BFGS);}
        SECTION("Beale")      {dlibTest(Beale, mini::algorithm::BFGS);}
        SECTION("decay2d")    {dlibTest(decay2d, mini::algorithm::BFGS);}
    }

    SECTION("dlib_global") {
        SECTION("decay") {dlibTest(decay1d, mini::algorithm::DLIB_GLOBAL);}
        SECTION("sqrt2") {dlibTest(sqrt2, mini::algorithm::DLIB_GLOBAL);}
        SECTION("euler") {dlibTest(euler, mini::algorithm::DLIB_GLOBAL);}
    }
}
#endif

TEST_CASE("Minimizer: levenberg_marquardt") {
    auto LMTest = [] (const TestFunction& test) {
        mini::LevenbergMarquardt mini(test.function, test.get_parameters());
        auto res = mini.minimize();
        CHECK(res.status == 0);
        for (int i = 0; i < static_cast<int>(test.min.size()); ++i) {
            CHECK_THAT(res.get_parameter(i).value, Catch::Matchers::WithinAbs(test.min[i], 1e-6));
        }
    };

    SECTION("decay")      {LMTest(decay1d);}
    SECTION("sqrt2")      {LMTest(sqrt2);}
    SECTION("euler")      {LMTest(euler);}
    SECTION("Rosenbrock") {LMTest(Rosenbrock);}
    SECTION("Beale")      {LMTest(Beale);}
    SECTION("decay2d")    {LMTest(decay2d);}

    SECTION("Rosenbrock from the classic start") {
        mini::LevenbergMarquardt mini(Rosenbrock.function, {mini::Parameter("a", -1.2, {-2, 2}), mini::Parameter("b", 1, {-2, 2})});
        auto res = mini.minimize();
        CHECK(res.status == 0);
        CHECK_THAT(res.get_parameter("a").value, Catch::Matchers::WithinAbs(1, 1e-6));
        CHECK_THAT(res.get_parameter("b").value, Catch::Matchers::WithinAbs(1, 1e-6));
    }

    SECTION("exponential decay with offset") {
        // noise-free y = 3 exp(-0.7 x) + 0.5, so the true parameters must be recovered exactly
        std::vector<double> x, y;
        for (int i = 0; i < 50; ++i) {x.push_back(0.1*i); y.push_back(3*std::exp(-0.7*x.back()) + 0.5);}
        mini::Minimizer::residual_function r = [&] (const std::vector<double>& p) {
            std::vector<double> res(x.size());
            for (int i = 0; i < static_cast<int>(x.size()); ++i) {res[i] = y[i] - (p[0]*std::exp(-p[1]*x[i]) + p[2]);}
            return res;
        };
        mini::LevenbergMarquardt mini(r, {mini::Parameter("A", 1, {0, 10}), mini::Parameter("k", 0.1, {0, 5}), mini::Parameter("c", 0, {-5, 5})});
        auto res = mini.minimize();
        CHECK(res.status == 0);
        CHECK_THAT(res.get_parameter("A").value, Catch::Matchers::WithinAbs(3, 1e-6));
        CHECK_THAT(res.get_parameter("k").value, Catch::Matchers::WithinAbs(0.7, 1e-6));
        CHECK_THAT(res.get_parameter("c").value, Catch::Matchers::WithinAbs(0.5, 1e-6));
    }

    SECTION("minimum outside the bounds") {
        // unconstrained minimum near (3, -1); the box pins the first parameter to its upper bound, while the second is free
        mini::Minimizer::residual_function r = [] (const std::vector<double>& p) {return std::vector{p[0] - 3, p[1] + 1, 0.1*p[0]*p[1]};};
        mini::LevenbergMarquardt mini(r, {mini::Parameter("a", 0.5, {0, 2}), mini::Parameter("b", 0, {-5, 5})});
        auto res = mini.minimize();
        CHECK(res.status == 0);
        CHECK_THAT(res.get_parameter("a").value, Catch::Matchers::WithinAbs(2, 1e-9));
        CHECK_THAT(res.get_parameter("b").value, Catch::Matchers::WithinAbs(-1/1.04, 1e-6)); // minimizes (b+1)^2 + (0.2b)^2
    }
}

TEST_CASE("Minimizer: create_minimizer") {
    std::vector<mini::algorithm> algorithms = {mini::algorithm::GOLDEN, mini::algorithm::LEVENBERG_MARQUARDT};
    #ifdef DLIB_AVAILABLE
        algorithms.push_back(mini::algorithm::BFGS);
    #endif

    for (auto t : algorithms) {
        auto mini = mini::create_minimizer(t, euler.function, {"a", euler.min[0] - 1, euler.bounds[0]});
        auto res = mini->minimize();
        CHECK_THAT(res.get_parameter("a").value, Catch::Matchers::WithinAbs(euler.min[0], mini->tol));
    }
}
