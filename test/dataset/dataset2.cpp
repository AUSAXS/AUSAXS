#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <settings/All.h>
#include <dataset/Dataset.h>
#include <dataset/Dataset2D.h>
#include <data/Protein.h>
#include <data/Body.h>
#include <data/Water.h>
#include <data/Atom.h>
#include <em/ImageStack.h>
#include <plots/all.h>

#include <iostream>

using std::vector;

TEST_CASE("plots") {
    std::vector<double> x = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    std::vector<double> y = {-6, -4, -1, 2, 1, 3, 6, 7, 9};
    SimpleDataset data(x, y);

    plots::PlotDataset::quick_plot(data, "test.plot");
}

TEST_CASE("rebin", "[manual]") {
    settings::general::verbose = false;
    SimpleDataset data("data/SHOC2/SHOC2.dat");
    SimpleDataset data_unbinned = data;
    data.rebin();
    data.save("temp/dataset/rebin/rebinned.dat");

    data.add_plot_options(style::draw::errors, {{"color", style::color::orange}, {"logx", true}, {"logy", true}});
    data_unbinned.add_plot_options(style::draw::errors, {{"color", style::color::black}, {"logx", true}, {"logy", true}});

    plots::PlotDataset plot(data_unbinned);
    plot.plot(data);
    plot.save("temp/dataset/rebin/both.png");

    plots::PlotDataset::quick_plot(data_unbinned, "temp/dataset/rebin/original.png");
    plots::PlotDataset::quick_plot(data, "temp/dataset/rebin/rebinned.png");
}

TEST_CASE("sim_err", "[manual]") {
    settings::general::verbose = false;
    Dataset2D data1("test/files/2epe.dat");
    Dataset2D data2 = data1;

    data2.simulate_errors();

    data1.add_plot_options(style::draw::points, {{"color", style::color::black}, {"lw", 2}});
    data2.add_plot_options(style::draw::points, {{"color", style::color::orange}, {"lw", 2}});

    plots::PlotIntensity plot(data1);
    plot.plot(data2);
    plot.save("temp/dataset/compare_errors.png");
}

TEST_CASE("sim_noise", "[manual]") {
    SimpleDataset data = SimpleDataset::generate_random_data(10000, 0, 1);
    SimpleDataset data2(data);
    data.simulate_noise();

    hist::Histogram hist;
    for (unsigned int i = 0; i < data.size(); i++) {
        double diff = (data2.y(i) - data.y(i))/data2.yerr(i);
        hist.p.push_back(diff);
    }

    hist.generate_axis();
    plots::PlotHistogram::quick_plot(hist, "temp/dataset/gaussian_noise.png");
}

TEST_CASE("reduce") {
    vector<double> x = {1, 2, 3, 4, 5};
    vector<double> y = {10, 20, 30, 40, 50};
    Dataset2D data(x, y, "i", "j");

    SECTION("reduce") {
        data.reduce(2);
        CHECK(data.size() < x.size());
    }

    SECTION("logreduce") {
        x = {1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9};
        y = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        data = Dataset2D(x, y);
        data.reduce(5, true);
        CHECK(data.size() == 5);
        CHECK(data.x() == std::vector<double>({1e0, 1e2, 1e4, 1e6, 1e8}));
    }
}

TEST_CASE("ranges") {
    std::vector<double> x = {1, 2, 3, 4, 5};
    std::vector<double> y = {10, 20, 30, 40, 50};
    SimpleDataset data(x, y);

    SECTION("limit") {
        SECTION("x") {
            data.limit_x(Limit(2, 3));
            CHECK(data.x() == std::vector<double>{2, 3});
            CHECK(data.y() == std::vector<double>{20, 30});
        }

        SECTION("y") {
            data.limit_y(Limit(20, 40));
            CHECK(data.x() == std::vector<double>{2, 3, 4});
            CHECK(data.y() == std::vector<double>{20, 30, 40});
        }
    }

    SECTION("spans") {
        SECTION("x") {
            auto span = data.span_x();
            CHECK(span == Limit(1, 5));
            CHECK(span == data.get_xlimits());
        }

        SECTION("y") {
            auto span = data.span_y();
            CHECK(span == Limit(10, 50));
            CHECK(span == data.get_ylimits());
        }

        SECTION("positive y") {
            y = {-6, -2, 1, 5, 8};
            data = SimpleDataset(x, y);
            auto span = data.span_y_positive();
            CHECK(span == Limit(1, 8));

            data = SimpleDataset();
            span = data.span_y_positive();
            CHECK(span == Limit(0, 0));
        }
    }
}

TEST_CASE("point2d") {
    Point2D p0(0, 4);
    Point2D p1(1, 2);
    Point2D p2(3, 4);
    Point2D p3(5, 6);
    Point2D p4(7, 8);
    Point2D p5(9, 10);

    Dataset2D data;
    data.push_back(p0);
    data.push_back(p1);
    data.push_back(p2);
    data.push_back(p3);
    data.push_back(p4);
    data.push_back(p5);

    CHECK(data.size() == 6);
    CHECK(data.x() == std::vector<double>({0, 1, 3, 5, 7, 9}));
    CHECK(data.y() == std::vector<double>({4, 2, 4, 6, 8, 10}));

    CHECK(data.get_point(0) == p0);
    CHECK(data.get_point(1) == p1);
    CHECK(data.get_point(2) == p2);

    CHECK(data.find_minimum() == p1);
}

TEST_CASE("reduceplot", "[manual]") {
    Protein protein("test/files/2epe.pdb");
    auto h = protein.get_histogram();

    plots::PlotIntensity plot(h);
    SimpleDataset data = h.calc_debye_scattering_intensity();
    data.reduce(20);
    plot.plot(data);
    plot.save("reduce_test.png");
}

// TEST_CASE("moving_average_plot", "[manual]") {
//     vector<double> x, y;
//     for (double xx = 0; xx < 2*M_PI; xx += 0.05) {
//         x.push_back(xx);
//         y.push_back(sin(xx));
//     }
//     SimpleDataset data(x, y, vector<double>(x.size(), 1));
//     data.add_plot_options("points");
//     plots::PlotDataset plot(data);

//     data = data.rolling_average(5);
//     data.interpolate(5);
//     data.add_plot_options(style::draw::line, {{"color", style::color::red}});
//     plot.plot(data);
//     plot.save("figures/test/dataset/moving_average.png");
// }