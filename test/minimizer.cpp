#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <dlib/optimization.h>
#include <dlib/global_optimization.h>

#include <mini/all.h>
#include <plots/all.h>

using std::vector;

struct TestFunction {
    TestFunction(std::function<double(std::vector<double>)> function, std::vector<Limit> bounds, std::vector<double> min) : function(function), bounds(bounds), min(min) {}
    TestFunction(std::function<double(std::vector<double>)> function, Limit bounds, double min) : TestFunction(function, vector{bounds}, vector{min}) {}

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
    return (1 - 8*x1 + 7*x1*x1 - 7*std::pow(x1, 3)/3 + std::pow(x1, 4)/4)*x2*x2*std::exp(-x1);}, 
    {Limit(0, 5), Limit(0, 5)}, 
    {4, 2}
);
TestFunction RosenbrockModified([] (std::vector<double> pars) {
    double x1 = pars[0], x2 = pars[1]; 
    return 74 + 100*std::pow(x2 - x1*x1, 2) + std::pow(1-x1, 2) - 400*std::exp(-(std::pow(x1+1, 2)+std::pow(x2+1, 2))*10);}, 
    {Limit(-1, 0.4), Limit(-1, 1)}, 
    {-0.9, -0.95}
);

TEST_CASE("1d_landscape", "[minimizer],[manual]") {
    mini::Golden mini(problem04.function, {"x1", problem04.bounds[0]});
    auto res = mini.minimize();

    SimpleDataset evaluations = mini.get_evaluated_points().as_dataset();
    SimpleDataset landscape = mini.landscape().as_dataset();
    landscape.add_plot_options(style::draw::line, {{"color", style::color::black}});
    evaluations.add_plot_options(style::draw::points, {{"color", style::color::orange}});

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
        CHECK_THAT(res.get_parameter("a").value, Catch::Matchers::WithinAbs(test.min[0], mini.tol));
    };

    SECTION("problem04") {GoldenTest(problem04);}
    SECTION("problem13") {GoldenTest(problem13);}
    SECTION("problem18") {GoldenTest(problem18);}
}

TEST_CASE("scan_minimizer", "[minimizer]") {
    auto ScanTest1D = [] (const TestFunction& test) {
        mini::Scan mini(test.function, {"a", test.bounds[0]});
        mini.set_evals(1000);
        auto res = mini.minimize();
        CHECK_THAT(res.get_parameter("a").value, Catch::Matchers::WithinAbs(test.min[0], mini.tol));
    };

    auto ScanTest1DRough = [] (const TestFunction& test) {
        mini::Scan mini(test.function, {"a", test.bounds[0]});
        mini.set_evals(10);
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

TEST_CASE("minimum_explorer", "[minimizer],[manual]") {
    auto ExplorerTest1D = [] (const TestFunction& test) {
        mini::dlibMinimizer mini1(test.function, {"a", test.bounds[0]});
        auto res = mini1.minimize();

        mini::Parameter p = res.get_parameter("a");
        mini::MinimumExplorer mini2(test.function, p, 100);
        res = mini2.minimize();

        SimpleDataset data1 = mini1.get_evaluated_points().as_dataset();
        data1.add_plot_options(style::draw::points, {{"xlabel", "x"}, {"ylabel", "f(x)"}, {"color", style::color::blue}});
        plots::PlotDataset plot(data1);

        SimpleDataset data2 = mini2.get_evaluated_points().as_dataset();
        data2.add_plot_options(style::draw::points, {{"xlabel", "x"}, {"ylabel", "f(x)"}, {"color", style::color::orange}});
        plots::PlotDataset::quick_plot(data2, "figures/test/minimizer/explorer_test_single.pdf");
        plot.plot(data2);

        mini::Golden mini3(test.function, {"a", test.bounds[0]});
        SimpleDataset line = mini3.landscape(1000).as_dataset();
        plot.plot(line);

        plot.save("figures/test/minimizer/explorer_test.pdf");

        CHECK_THAT(res.get_parameter("a").value, Catch::Matchers::WithinAbs(test.min[0], mini1.tol));
    };

    // test with a fine grid
    SECTION("problem04") {ExplorerTest1D(problem04);}
    // SECTION("problem13") {ExplorerTest1D(problem13);}
    // SECTION("problem18") {ExplorerTest1D(problem18);}
}

typedef dlib::matrix<double,0,1> column_vector;
TEST_CASE("dlib", "[minimizer]") {
    auto dlibTest1D = [] (const TestFunction& test, auto strat) {
        auto wrapper = [&test] (const column_vector& x) {
            return test.function({x(0)});
        };

        column_vector starting_point = {test.get_center()[0]};
        dlib::find_min_box_constrained(strat, dlib::objective_delta_stop_strategy(1e-7), wrapper, dlib::derivative(wrapper), starting_point, test.bounds[0].min, test.bounds[0].max);
        CHECK_THAT(starting_point(0), Catch::Matchers::WithinAbs(test.min[0], 1e-3));
    };

    auto dlibTest2D = [] (const TestFunction& test, auto strat) {
        auto wrapper = [&test] (const column_vector& x) {
            std::vector<double> v = {x(0), x(1)};
            return test.function(v);
        };

        column_vector starting_point = {test.get_center()[0], test.get_center()[1]};
        column_vector bmin = {test.bounds[0].min, test.bounds[1].min};
        column_vector bmax = {test.bounds[0].max, test.bounds[1].max};
        dlib::find_min_box_constrained(strat, dlib::objective_delta_stop_strategy(1e-7), wrapper, dlib::derivative(wrapper), starting_point, bmin, bmax);

        CHECK_THAT(starting_point(0), Catch::Matchers::WithinAbs(test.min[0], 1e-3));
        CHECK_THAT(starting_point(1), Catch::Matchers::WithinAbs(test.min[1], 1e-3));
    };

    SECTION("bfgs") {
        SECTION("problem04") {dlibTest1D(problem04, dlib::bfgs_search_strategy());}
        SECTION("problem13") {dlibTest1D(problem13, dlib::bfgs_search_strategy());}
        SECTION("problem18") {dlibTest1D(problem18, dlib::bfgs_search_strategy());}

        SECTION("Decanomial") {dlibTest2D(Decanomial, dlib::bfgs_search_strategy());}
        SECTION("Hosaki") {dlibTest2D(Hosaki, dlib::bfgs_search_strategy());}
        SECTION("RosenbrockModified") {dlibTest2D(RosenbrockModified, dlib::bfgs_search_strategy());}
    }
}

TEST_CASE("dlib_interface", "[minimizer]") {
    auto dlibTest1D = [] (const TestFunction& test) {
        mini::dlibMinimizer mini(test.function, {"a", test.bounds[0]});
        auto res = mini.minimize();
        CHECK_THAT(res.get_parameter("a").value, Catch::Matchers::WithinAbs(test.min[0], mini.tol));
    };
    SECTION("problem04") {dlibTest1D(problem04);}
    SECTION("problem13") {dlibTest1D(problem13);}
    SECTION("problem18") {dlibTest1D(problem18);}
}

TEST_CASE("create_minimizer", "[minimizer]") {
    SECTION("dlib") {
        auto dlib = mini::create_minimizer(mini::type::BFGS, problem04.function, {"a", problem04.bounds[0]});
        auto res = dlib->minimize();
        CHECK_THAT(res.get_parameter("a").value, Catch::Matchers::WithinAbs(problem04.min[0], dlib->tol));
    }

    SECTION("golden") {
        auto golden = mini::create_minimizer(mini::type::GOLDEN, problem04.function, {"a", problem04.bounds[0]});
        auto res = golden->minimize();
        CHECK_THAT(res.get_parameter("a").value, Catch::Matchers::WithinAbs(problem04.min[0], golden->tol));
    }
}