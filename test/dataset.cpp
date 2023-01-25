#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <utility/Settings.h>
#include <utility/Dataset.h>
#include <em/ImageStack.h>
#include <plots/all.h>

#include <iostream>

using std::vector;

TEST_CASE("dataset_basics", "[dataset]") {
    std::vector<double> xd = {   1,   2,   3,   4,   5,   6,   7,   8,   9};
    std::vector<double> yd = {  -6,  -4,  -1,   2,   1,   3,   6,   7,   9};
    std::vector<double> xed = {0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5};
    std::vector<double> yed = { -7,  -5,  -2,   1,   0,   2,   5,   6,   8};
    Dataset2D data(xd, yd, xed, yed);

    SECTION("columns") {
        auto x = data.x();
        auto y = data.y();
        auto yerr = data.yerr();
        auto xerr = data.xerr();

        REQUIRE(x.size() == 9);
        REQUIRE(y.size() == 9);
        REQUIRE(yerr.size() == 9);
        REQUIRE(xerr.size() == 9);
        CHECK(x == Vector(xd));
        CHECK(y == Vector(yd));
        CHECK(yerr == Vector(yed));
        CHECK(xerr == Vector(xed));

        // columns are mutable
        std::transform(x.begin(), x.end(), yerr.begin(), [](double x) {return x;});
        std::transform(y.begin(), y.end(), yerr.begin(), y.begin(), std::multiplies<>());
        CHECK(y == Vector{-6, -8, -3, 8, 5, 18, 42, 56, 81});
    }

    SECTION("column by name") {
        std::vector<std::string> names = {"x", "y", "yerr", "xerr"};
        data.set_col_names(names);
        auto x = data.col("x");
        auto y = data.col("y");
        auto yerr = data.col("yerr");
        auto xerr = data.col("xerr");

        REQUIRE(data.get_col_names() == names);
        CHECK(data.get_col_names(0) == "x");
        CHECK(data.get_col_names(1) == "y");
        CHECK(data.get_col_names(2) == "yerr");
        CHECK(data.get_col_names(3) == "xerr");

        CHECK(x == Vector(xd));
        CHECK(y == Vector(yd));
        CHECK(yerr == Vector(yed));
        CHECK(xerr == Vector(xed));
    }

    SECTION("indexing") {
        CHECK(data.x(0) == 1);
        CHECK(data.x(1) == 2);
        CHECK(data.x(2) == 3);

        CHECK(data.y(0) == -6);
        CHECK(data.y(1) == -4);
        CHECK(data.y(2) == -1);

        CHECK(data.xerr(0) == 0.5);
        CHECK(data.xerr(1) == 1.5);
        CHECK(data.xerr(2) == 2.5);

        CHECK(data.yerr(0) == -7);
        CHECK(data.yerr(1) == -5);
        CHECK(data.yerr(2) == -2);

        const auto& cdata = data;
        CHECK(cdata.x() == data.x());
        CHECK(cdata.x(0) == 1);
        CHECK(cdata.x(1) == 2);
        CHECK(cdata.x(2) == 3);

        CHECK(cdata.y() == data.y());
        CHECK(cdata.y(0) == -6);
        CHECK(cdata.y(1) == -4);
        CHECK(cdata.y(2) == -1);

        CHECK(cdata.xerr() == data.xerr());
        CHECK(cdata.xerr(0) == 0.5);
        CHECK(cdata.xerr(1) == 1.5);
        CHECK(cdata.xerr(2) == 2.5);

        CHECK(cdata.yerr() == data.yerr());
        CHECK(cdata.yerr(0) == -7);
        CHECK(cdata.yerr(1) == -5);
        CHECK(cdata.yerr(2) == -2);
    }

    SECTION("constructors") {
        Dataset2D data2(data);
        CHECK(data2.x() == data.x());
        CHECK(data2.y() == data.y());
        CHECK(data2.xerr() == data.xerr());
        CHECK(data2.yerr() == data.yerr());

        Dataset2D data3(std::move(data));
        CHECK(data3.x() == data2.x());
        CHECK(data3.y() == data2.y());
        CHECK(data3.xerr() == data2.xerr());
        CHECK(data3.yerr() == data2.yerr());

        Dataset data4(vector<vector<double>>{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
        CHECK(data4.x() == Vector{1, 2, 3});
        CHECK(data4.y() == Vector{4, 5, 6});
        CHECK(data4.col(2) == Vector{7, 8, 9});
        SimpleDataset data5(data4);
        CHECK(data5.x() == Vector{1, 2, 3});
        CHECK(data5.y() == Vector{4, 5, 6});
        CHECK(data5.yerr() == Vector{7, 8, 9});
        CHECK(data5.N == 3);
        CHECK(data5.M == 3);

        data4 = Dataset(vector<vector<double>>{{1, 3, 5, 7, 9}, {2, 4, 6, 8, 10}});
        CHECK(data4.x() == Vector{1, 3, 5, 7, 9});
        CHECK(data4.y() == Vector{2, 4, 6, 8, 10});
        SimpleDataset data6(data4);
        CHECK(data6.x() == Vector{1, 3, 5, 7, 9});
        CHECK(data6.y() == Vector{2, 4, 6, 8, 10});
        CHECK(data6.N == 5);
        CHECK(data6.M == 3);
    }
}

TEST_CASE("dataset_pushback", "[dataset]") {
    SECTION("Dataset2D") {
        Dataset2D data;
        data.push_back(1, 2, 3, 4);
        data.push_back(5, 6, 7, 8);
        data.push_back(9, 10, 11, 12);
        data.push_back(Point2D(13, 14, 15, 16));        

        CHECK(data.x() == Vector{1, 5, 9, 13});
        CHECK(data.y() == Vector{2, 6, 10, 14});
        CHECK(data.xerr() == Vector{3, 7, 11, 15});
        CHECK(data.yerr() == Vector{4, 8, 12, 16});
    }

    SECTION("SimpleDataset") {
        SimpleDataset data;
        data.push_back(1, 2);
        data.push_back(3, 4);
        data.push_back(5, 6);
        data.push_back(Point2D(7, 8));

        CHECK(data.x() == Vector{1, 3, 5, 7});
        CHECK(data.y() == Vector{2, 4, 6, 8});
    }
}

TEST_CASE("dataset_ylimits", "[dataset]") {
    std::vector<double> x = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    std::vector<double> y = {-6, -4, -1, 2, 1, 3, 6, 7, 9};
    SimpleDataset data(x, y);

    Limit limit(0.5, 5);
    data.limit_y(limit);
    for (unsigned int i = 0; i < data.size(); i++) {
        CHECK(limit.min <= data.y(i));
        CHECK(data.y(i) <= limit.max);
    }
}

TEST_CASE("dataset_xlimits", "[dataset]") {
    std::vector<double> x = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    std::vector<double> y = {-6, -4, -1, 2, 1, 3, 6, 7, 9};
    SimpleDataset data(x, y);

    Limit limit(0.5, 5);
    data.limit_x(limit);
    for (unsigned int i = 0; i < data.size(); i++) {
        CHECK(limit.min <= data.x(i));
        CHECK(data.x(i) <= limit.max);
    }
}

TEST_CASE("dataset_plots", "[dataset]") {
    std::vector<double> x = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    std::vector<double> y = {-6, -4, -1, 2, 1, 3, 6, 7, 9};
    SimpleDataset data(x, y);

    plots::PlotDataset::quick_plot(data, "test.plot");
}

TEST_CASE("dataset_rebin", "[dataset],[files],[manual]") {
    setting::general::verbose = false;
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

TEST_CASE("dataset_sim_err", "[dataset],[files],[manual]") {
    setting::general::verbose = false;
    Dataset2D data1("test/files/2epe.dat");
    Dataset2D data2 = data1;

    data2.simulate_errors();

    data1.add_plot_options(style::draw::points, {{"color", style::color::black}, {"lw", 2}});
    data2.add_plot_options(style::draw::points, {{"color", style::color::orange}, {"lw", 2}});

    plots::PlotIntensity plot(data1);
    plot.plot(data2);
    plot.save("temp/dataset/compare_errors.png");
}

TEST_CASE("dataset_sim_noise", "[dataset],[manual]") {
    SimpleDataset data = SimpleDataset::generate_random_data(10000);
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

TEST_CASE("dataset_io", "[dataset],[files]") {
    setting::general::verbose = false;
    SECTION("lysozyme") {
        Dataset2D data("test/files/2epe.dat");
        data.save("temp/dataset/2epe.dat");
        Dataset2D data2("temp/dataset/2epe.dat");
        REQUIRE(data.size() == data2.size());
        REQUIRE(data.x().size() == data2.x().size());
        REQUIRE(data.y().size() == data2.y().size());
        REQUIRE(data.yerr().size() == data2.yerr().size());

        for (unsigned int i = 0; i < data.size(); i++) {
            CHECK_THAT(data.x(i), Catch::Matchers::WithinRel(data2.x(i), 1e-3));
            CHECK_THAT(data.y(i), Catch::Matchers::WithinRel(data2.y(i), 1e-3));
            CHECK_THAT(data.yerr(i), Catch::Matchers::WithinRel(data2.yerr(i), 1e-3));
        }
    }
}

TEST_CASE("dataset_read", "[dataset],[files]") {
    setting::general::verbose = false;
    SECTION("actual data") {
        Dataset2D data("test/files/2epe.dat");
        auto x = data.x();
        auto y = data.y();
        auto yerr = data.yerr();

        vector<double> validate_x = {9.81300045E-03, 1.06309997E-02, 1.14489999E-02, 1.22659998E-02, 1.30840000E-02, 1.39020002E-02, 1.47200003E-02, 1.55379996E-02, 1.63550004E-02, 1.71729997E-02};
        vector<double> validate_y = {6.67934353E-03, 7.27293547E-03, 8.74083303E-03, 9.22449585E-03, 9.13867634E-03, 9.21153929E-03, 9.37998667E-03, 8.67372658E-03, 9.23649967E-03, 9.22480784E-03};
        vector<double> validate_yerr = {1.33646582E-03, 1.01892441E-03, 8.62116576E-04, 7.71059655E-04, 6.87870081E-04, 6.30189374E-04, 4.98525158E-04, 4.69041377E-04, 4.46073769E-04, 4.26004088E-04};

        REQUIRE(x.size() == 104);
        REQUIRE(y.size() == 104);
        REQUIRE(yerr.size() == 104);
        for (unsigned int i = 0; i < validate_x.size(); i++) {
            CHECK_THAT(x[i], Catch::Matchers::WithinRel(validate_x[i]));
            CHECK_THAT(y[i], Catch::Matchers::WithinRel(validate_y[i]));
            CHECK_THAT(yerr[i], Catch::Matchers::WithinRel(validate_yerr[i]));
        }
    }
}

TEST_CASE("dataset_scale", "[dataset]") {
    vector<double> x = {1, 2, 3, 4, 5};
    vector<double> y = {10, 20, 30, 40, 50};
    vector<double> yerr = {1, 2, 3, 4, 5};
    vector<double> xerr = {0.1, 0.2, 0.3, 0.4, 0.5};
    Dataset2D data(x, y, xerr, yerr);

    SECTION("scale") {
        data.scale_y(2);
        CHECK(data.y().to_vector() == vector<double>({20, 40, 60, 80, 100}));

        data.scale_y(0.2);
        CHECK(data.y().to_vector() == vector<double>({4, 8, 12, 16, 20}));

        data.scale_y(-0.5);
        CHECK(data.y().to_vector() == vector<double>({-2, -4, -6, -8, -10}));
    }

    SECTION("scale errors") {
        data.scale_errors(2);
        CHECK(data.yerr().to_vector() == vector<double>({2, 4, 6, 8, 10}));
        CHECK(data.xerr().to_vector() == vector<double>({0.2, 0.4, 0.6, 0.8, 1.0}));

        data.scale_errors(0.5);
        CHECK(data.yerr().to_vector() == vector<double>({1, 2, 3, 4, 5}));
        CHECK(data.xerr().to_vector() == vector<double>({0.1, 0.2, 0.3, 0.4, 0.5}));

        data.scale_errors(-0.5);
        CHECK(data.yerr().to_vector() == vector<double>({-0.5, -1, -1.5, -2, -2.5}));
        CHECK(data.xerr().to_vector() == vector<double>({-0.05, -0.1, -0.15, -0.2, -0.25}));
    }

    SECTION("normalize") {
        data.normalize(1);
        CHECK(data.y().to_vector() == vector<double>({1, 2, 3, 4, 5}));

        data.normalize(-5);
        CHECK(data.y().to_vector() == vector<double>({-5, -10, -15, -20, -25}));

        data.normalize(10);
        CHECK(data.y().to_vector() == vector<double>({10, 20, 30, 40, 50}));
    }
}

TEST_CASE("dataset_reduce", "[dataset]") {
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
        CHECK(data.x().to_vector() == vector<double>({1e0, 1e2, 1e4, 1e6, 1e8}));
    }
}

TEST_CASE("dataset_ranges", "[dataset]") {
    vector<double> x = {1, 2, 3, 4, 5};
    vector<double> y = {10, 20, 30, 40, 50};
    SimpleDataset data(x, y);

    SECTION("limit") {
        SECTION("x") {
            data.limit_x(Limit(2, 3));
            CHECK(data.x().to_vector() == vector<double>{2, 3});
            CHECK(data.y().to_vector() == vector<double>{20, 30});
        }

        SECTION("y") {
            data.limit_y(Limit(20, 40));
            CHECK(data.x().to_vector() == vector<double>{2, 3, 4});
            CHECK(data.y().to_vector() == vector<double>{20, 30, 40});
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

TEST_CASE("dataset_point2d", "[dataset]") {
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
    CHECK(data.x().to_vector() == vector<double>({0, 1, 3, 5, 7, 9}));
    CHECK(data.y().to_vector() == vector<double>({4, 2, 4, 6, 8, 10}));

    CHECK(data.get_point(0) == p0);
    CHECK(data.get_point(1) == p1);
    CHECK(data.get_point(2) == p2);

    CHECK(data.find_minimum() == p1);
}

TEST_CASE("dataset_remove_duplicates", "[dataset]") {
    vector<double> x = {1, 2, 3, 4, 5, 5, 5, 6, 7, 8, 9};
    vector<double> y = {10, 20, 30, 40, 50, 50, 50, 60, 70, 80, 90};
    Dataset2D data(x, y);

    data.remove_consecutive_duplicates();
    CHECK(data.size() == 9);
    CHECK(data.x().to_vector() == vector<double>({1, 2, 3, 4, 5, 6, 7, 8, 9}));
    CHECK(data.y().to_vector() == vector<double>({10, 20, 30, 40, 50, 60, 70, 80, 90}));
}

TEST_CASE("dataset_sort", "[dataset]") {
    vector<double> x = {1, 0, 4, 3, 2};
    vector<double> y = {10, 0, 40, 30, 20};
    Dataset2D data(x, y);

    data.sort_x();
    CHECK(data.x().to_vector() == vector<double>({0, 1, 2, 3, 4}));
    CHECK(data.y().to_vector() == vector<double>({0, 10, 20, 30, 40}));
}

TEST_CASE("dataset_reduceplot", "[dataset],[manual]") {
    Protein protein("test/files/2epe.pdb");
    auto h = protein.get_histogram();

    plots::PlotIntensity plot(h);
    SimpleDataset data = h.calc_debye_scattering_intensity();
    data.reduce(20);
    plot.plot(data);
    plot.save("reduce_test.png");
}

TEST_CASE("dataset_stats", "[dataset]") {
    Dataset2D data1(
        vector<double>{1, 1, 1, 1}, 
        vector<double>{10, 3, 5, 6}, 
        vector<double>{1, 2, 3, 4}
    );

    Dataset2D data2(
        vector<double>{1, 1, 1, 1, 1, 1}, 
        vector<double>{12, 14, 15, 15, 14, 17}, 
        vector<double>{1, 0.5, 0.1, 2, 0.9, 3}
    );

    Dataset2D data3(
        vector<double>{1, 1, 1, 1, 1, 1, 1, 1, 1}, 
        vector<double>{54, 66, 78, 80, 82, 84, 84, 90, 93}, 
        vector<double>{1, 0.2, 2, 0.5, 0.9, 4, 1, 0.1, 0.4}
    );

    SECTION("weighted_mean_error") {
        CHECK_THAT(data1.weighted_mean_error(), Catch::Matchers::WithinAbs(0.838116, 1e-6));
        CHECK_THAT(data2.weighted_mean_error(), Catch::Matchers::WithinAbs(0.0968561, 1e-6));
        CHECK_THAT(data3.weighted_mean_error(), Catch::Matchers::WithinAbs(0.0848808, 1e-6));
    }

    SECTION("weighted_mean") {
        CHECK_THAT(data1.weighted_mean(), Catch::Matchers::WithinAbs(8.204903, 1e-3));
        CHECK_THAT(data2.weighted_mean(), Catch::Matchers::WithinAbs(14.924834, 1e-3));
        CHECK_THAT(data3.weighted_mean(), Catch::Matchers::WithinAbs(85.125966, 1e-3));       
    }

    SECTION("mean") {
        CHECK_THAT(data1.mean(), Catch::Matchers::WithinAbs(6, 1e-3));
        CHECK_THAT(data2.mean(), Catch::Matchers::WithinAbs(14.5, 1e-3));
        CHECK_THAT(data3.mean(), Catch::Matchers::WithinAbs(79, 1e-3));
    }

    SECTION("std") {
        CHECK_THAT(data1.std(), Catch::Matchers::WithinAbs(2.943920, 1e-3));
        CHECK_THAT(data2.std(), Catch::Matchers::WithinAbs(1.643167, 1e-3));
        CHECK_THAT(data3.std(), Catch::Matchers::WithinAbs(12.103718, 1e-3));
    }
}

#include <math/MovingAverager.h>
TEST_CASE("dataset_moving_average", "[dataset]") {
    SECTION("simple") {
        SimpleDataset data(
            vector<double>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, 
            vector<double>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, 
            vector<double>{1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
        );
        SECTION("unweighted") {
            SECTION("3") {
                auto res = MovingAverage::average(data.y(), 3);
                REQUIRE(res.size() == 10);
                CHECK(res[0] == 1);
                CHECK(res[1] == 2);
                CHECK(res[2] == 3);
                CHECK(res[3] == 4);
                CHECK(res[4] == 5);
                CHECK(res[5] == 6);
                CHECK(res[6] == 7);
                CHECK(res[7] == 8);
                CHECK(res[8] == 9);
                CHECK(res[9] == 10);
            }

            SECTION("5") {
                auto res = MovingAverage::average(data.y(), 5);
                REQUIRE(res.size() == 10);
                CHECK(res[0] == 1);
                CHECK(res[1] == 2);
                CHECK(res[2] == 3);
                CHECK(res[3] == 4);
                CHECK(res[4] == 5);
                CHECK(res[5] == 6);
                CHECK(res[6] == 7);
                CHECK(res[7] == 8);
                CHECK(res[8] == 9);
                CHECK(res[9] == 10);
            }
        }

        SECTION("half_moving_average") {
            SECTION("3") {
                SimpleDataset res = data.rolling_average(3);
                REQUIRE(res.size() == 10);
                CHECK_THAT(res.x(0), Catch::Matchers::WithinAbs(1, 1e-6));
                CHECK_THAT(res.y(0), Catch::Matchers::WithinAbs(1, 1e-6));

                CHECK_THAT(res.x(1), Catch::Matchers::WithinAbs(2, 1e-6));
                CHECK_THAT(res.y(1), Catch::Matchers::WithinAbs((1./2 + 2 + 3./2)/2, 1e-6));

                CHECK_THAT(res.x(2), Catch::Matchers::WithinAbs(3, 1e-6));
                CHECK_THAT(res.y(2), Catch::Matchers::WithinAbs((2./2 + 3 + 4./2)/2, 1e-6));

                CHECK_THAT(res.x(3), Catch::Matchers::WithinAbs(4, 1e-6));
                CHECK_THAT(res.y(3), Catch::Matchers::WithinAbs((3./2 + 4 + 5./2)/2, 1e-6));

                CHECK_THAT(res.x(4), Catch::Matchers::WithinAbs(5, 1e-6));
                CHECK_THAT(res.y(4), Catch::Matchers::WithinAbs((4./2 + 5 + 6./2)/2, 1e-6));

                CHECK_THAT(res.x(5), Catch::Matchers::WithinAbs(6, 1e-6));
                CHECK_THAT(res.y(5), Catch::Matchers::WithinAbs((5./2 + 6 + 7./2)/2, 1e-6));

                CHECK_THAT(res.x(6), Catch::Matchers::WithinAbs(7, 1e-6));
                CHECK_THAT(res.y(6), Catch::Matchers::WithinAbs((6./2 + 7 + 8./2)/2, 1e-6));

                CHECK_THAT(res.x(7), Catch::Matchers::WithinAbs(8, 1e-6));
                CHECK_THAT(res.y(7), Catch::Matchers::WithinAbs((7./2 + 8 + 9./2)/2, 1e-6));

                CHECK_THAT(res.x(8), Catch::Matchers::WithinAbs(9, 1e-6));
                CHECK_THAT(res.y(8), Catch::Matchers::WithinAbs((8./2 + 9 + 10./2)/2, 1e-6));

                CHECK_THAT(res.x(9), Catch::Matchers::WithinAbs(10, 1e-6));
                CHECK_THAT(res.y(9), Catch::Matchers::WithinAbs(10, 1e-6));
            }

            SECTION("5") {
                SimpleDataset res = data.rolling_average(5);
                REQUIRE(res.size() == 10);
                CHECK_THAT(res.x(0), Catch::Matchers::WithinAbs(1, 1e-6));
                CHECK_THAT(res.y(0), Catch::Matchers::WithinAbs(1, 1e-6));

                CHECK_THAT(res.x(1), Catch::Matchers::WithinAbs(2, 1e-6));
                CHECK_THAT(res.y(1), Catch::Matchers::WithinAbs((1./2 + 2 + 3./2)/2, 1e-6));

                CHECK_THAT(res.x(2), Catch::Matchers::WithinAbs(3, 1e-6));
                CHECK_THAT(res.y(2), Catch::Matchers::WithinAbs((1./4 + 2./2 + 3 + 4./2 + 5./4)/2.5, 1e-6));

                CHECK_THAT(res.x(3), Catch::Matchers::WithinAbs(4, 1e-6));
                CHECK_THAT(res.y(3), Catch::Matchers::WithinAbs((2./4 + 3./2 + 4 + 5./2 + 6./4)/2.5, 1e-6));

                CHECK_THAT(res.x(4), Catch::Matchers::WithinAbs(5, 1e-6));
                CHECK_THAT(res.y(4), Catch::Matchers::WithinAbs((3./4 + 4./2 + 5 + 6./2 + 7./4)/2.5, 1e-6));

                CHECK_THAT(res.x(5), Catch::Matchers::WithinAbs(6, 1e-6));
                CHECK_THAT(res.y(5), Catch::Matchers::WithinAbs((4./4 + 5./2 + 6 + 7./2 + 8./4)/2.5, 1e-6));

                CHECK_THAT(res.x(6), Catch::Matchers::WithinAbs(7, 1e-6));
                CHECK_THAT(res.y(6), Catch::Matchers::WithinAbs((5./4 + 6./2 + 7 + 8./2 + 9./4)/2.5, 1e-6));

                CHECK_THAT(res.x(7), Catch::Matchers::WithinAbs(8, 1e-6));
                CHECK_THAT(res.y(7), Catch::Matchers::WithinAbs((6./4 + 7./2 + 8 + 9./2 + 10./4)/2.5, 1e-6));

                CHECK_THAT(res.x(8), Catch::Matchers::WithinAbs(9, 1e-6));
                CHECK_THAT(res.y(8), Catch::Matchers::WithinAbs((8./2 + 9 + 10./2)/2, 1e-6));

                CHECK_THAT(res.x(9), Catch::Matchers::WithinAbs(10, 1e-6));
                CHECK_THAT(res.y(9), Catch::Matchers::WithinAbs(10, 1e-6));
            }
        }
    }
}

TEST_CASE("dataset_interpolate", "[dataset]") {
    SECTION("simple") {
        SimpleDataset data(
            vector<double>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, 
            vector<double>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}
        );

        data.interpolate(1);
        REQUIRE(data.size() == 18);
        CHECK_THAT(data.x(0), Catch::Matchers::WithinAbs(1, 1e-6));
        CHECK_THAT(data.y(0), Catch::Matchers::WithinAbs(1, 1e-6));

        CHECK_THAT(data.x(1), Catch::Matchers::WithinAbs(1.5, 1e-6));
        CHECK_THAT(data.y(1), Catch::Matchers::WithinAbs(1.5, 1e-6));

        CHECK_THAT(data.x(2), Catch::Matchers::WithinAbs(2, 1e-6));
        CHECK_THAT(data.y(2), Catch::Matchers::WithinAbs(2, 1e-6));

        CHECK_THAT(data.x(3), Catch::Matchers::WithinAbs(2.5, 1e-6));
        CHECK_THAT(data.y(3), Catch::Matchers::WithinAbs(2.5, 1e-6));

        CHECK_THAT(data.x(4), Catch::Matchers::WithinAbs(3, 1e-6));
        CHECK_THAT(data.y(4), Catch::Matchers::WithinAbs(3, 1e-6));

        CHECK_THAT(data.x(5), Catch::Matchers::WithinAbs(3.5, 1e-6));
        CHECK_THAT(data.y(5), Catch::Matchers::WithinAbs(3.5, 1e-6));
    }

    SECTION("sine") {
        vector<double> x, y;
        for (double xx = 0; xx < 2*M_PI; xx += 0.05) {
            x.push_back(xx);
            y.push_back(sin(xx));
        }

        SimpleDataset data(x, y);
        data.interpolate(5);
        for (unsigned int i = 0; i < data.size(); i++) {
            CHECK_THAT(data.y(i), Catch::Matchers::WithinAbs(sin(data.x(i)), 1e-3));
        }
    }
}

TEST_CASE("dataset_moving_average_plot", "[dataset],[manual]") {
    vector<double> x, y;
    for (double xx = 0; xx < 2*M_PI; xx += 0.05) {
        x.push_back(xx);
        y.push_back(sin(xx));
    }
    SimpleDataset data(x, y, vector<double>(x.size(), 1));
    data.add_plot_options("points");
    plots::PlotDataset plot(data);

    data = data.rolling_average(5);
    data.interpolate(5);
    data.add_plot_options(style::draw::line, {{"color", style::color::red}});
    plot.plot(data);
    plot.save("figures/test/dataset/moving_average.png");
}

TEST_CASE("dataset_io_accuracy", "[dataset]") {
    SimpleDataset data;
    for (double x = 1.347e-01; x < 1.351e-01; x += 1e-6) {
        data.push_back(x, sin(x));
    }
    data.save("temp/dataset_io_accuracy.dat");

    SimpleDataset data2("temp/dataset_io_accuracy.dat");
    REQUIRE(data.size() == data2.size());
    for (unsigned int i = 0; i < data.size(); i++) {
        CHECK_THAT(data.x(i), Catch::Matchers::WithinAbs(data2.x(i), 1e-6));
        CHECK_THAT(data.y(i), Catch::Matchers::WithinAbs(data2.y(i), 1e-6));
    }
}