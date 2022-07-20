#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <utility/Dataset.h>
#include <plots/all.h>

#include <iostream>

TEST_CASE("debug", "[dataset],[disable]") {
    Multiset data("temp/multiset");
    std::cout << data[0].size() << std::endl;
    std::cout << data[1].size() << std::endl;
    std::cout << data[2].size() << std::endl;
    std::cout << data[3].size() << std::endl;
    std::cout << data[4].size() << std::endl;

    plots::PlotResolutionComparison plot_r(data);
    plot_r.save("figures/test/dataset/debug.pdf");
}

TEST_CASE("dataset_ylimits", "[dataset]") {
    Dataset2D data;
    data.x() = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    data.y() = {-10, -6, -4, -1, 2, 1, 3, 6, 7, 9};

    Limit limit(0.5, 5);
    data.limit_y(limit);
    auto y = data.y();
    for (unsigned int i = 0; i < data.size(); i++) {
        CHECK(limit.min <= y[i]);
        CHECK(y[i] <= limit.max);
    }
}

TEST_CASE("dataset_xlimits", "[dataset]") {
    Dataset2D data;
    data.x() = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    data.y() = {-10, -6, -4, -1, 2, 1, 3, 6, 7, 9};

    Limit limit(0.5, 5);
    data.limit_x(limit);
    auto x = data.x();
    for (unsigned int i = 0; i < data.size(); i++) {
        CHECK(limit.min <= x[i]);
        CHECK(x[i] <= limit.max);
    }
}

TEST_CASE("dataset_rebin", "[dataset],[files],[manual]") {
    SAXSDataset data("data/SHOC2/7sd0.dat");
    SAXSDataset data_unbinned = data;
    data.rebin();
    data.save("temp/dataset/rebin.dat");

    plots::PlotDataset plot(data_unbinned);
    plot.plot(data);
    plot.save("temp/dataset_rebin.pdf");
}

TEST_CASE("dataset_sim_err", "[dataset],[files],[manual]") {
    SAXSDataset data1("data/2epe.RSR");
    SAXSDataset data2 = data1;
    data2.simulate_errors();

    data1.add_plot_options("markers", {{"color", kBlack}, {"lw", 2}});
    data1.add_plot_options("markers", {{"color", kOrange+2}, {"lw", 2}});

    plots::PlotIntensity plot(data1);
    plot.plot_intensity(data2);
    plot.save("temp/compare_errors.pdf");
}

TEST_CASE("dataset_sim_noise", "[dataset],[manual]") {
    Dataset data1 = Dataset::generate_random_data(10000);
    Dataset data2(data1);
    data1.simulate_noise();

    hist::Histogram hist;
    auto y1 = data1.y();
    auto y2 = data2.y();
    auto yerr2 = data2.yerr();
    for (unsigned int i = 0; i < data1.size(); i++) {
        double diff = (y2[i] - y1[i])/yerr2[i];
        hist.p.push_back(diff);
    }

    hist.generate_axis();
    plots::PlotHistogram::quick_plot(hist, "temp/dataset/gaussian_noise.pdf");
}

TEST_CASE("dataset_is_logarithmic", "[dataset],[files]") {
    SECTION("lysozyme") {
        Dataset data("data/2epe.RSR");
        CHECK(data.is_logarithmic());
    }

    SECTION("A2M") {
        Dataset data("data/A2M_ma.RSR");
        CHECK(data.is_logarithmic());
    }

    SECTION("linear") {
        vector<double> x = {1, 2, 3, 4, 5, 6, 7, 8, 9};
        vector<double> y = {0, 1, 2, 3, 4, 5, 6, 7, 8};
        Dataset data(x, y);
        CHECK(!data.is_logarithmic());
    }
}

TEST_CASE("dataset_read", "[dataset],[files]") {
    SAXSDataset data("data/2epe.RSR");
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

TEST_CASE("dataset_normalize", "[dataset]") {
    vector<double> x = {1, 2, 3, 4, 5};
    vector<double> y = {10, 20, 30, 40, 50};
    Dataset data(x, y);

    SECTION("scale") {
        auto y = data.y();
        data.scale_y(2);
        CHECK(y[0] == 20);
        CHECK(y[1] == 40);
        CHECK(y[2] == 60);
        CHECK(y[3] == 80);
        CHECK(y[4] == 100);

        data.scale_y(0.2);
        CHECK(y[0] == 4);
        CHECK(y[1] == 8);
        CHECK(y[2] == 12);
        CHECK(y[3] == 16);
        CHECK(y[4] == 20);

        data.scale_y(-0.5);
        CHECK(y[0] == -2);
        CHECK(y[1] == -4);
        CHECK(y[2] == -6);
        CHECK(y[3] == -8);
        CHECK(y[4] == -10);
    }

    SECTION("normalize") {
        auto y = data.y();
        data.normalize(1);
        CHECK(y[0] == 1);
        CHECK(y[1] == 2);
        CHECK(y[2] == 3);
        CHECK(y[3] == 4);
        CHECK(y[4] == 5);

        data.normalize(-5);
        CHECK(y[0] == -5);
        CHECK(y[1] == -10);
        CHECK(y[2] == -15);
        CHECK(y[3] == -20);
        CHECK(y[4] == -25);

        data.normalize(10);
        CHECK(y[0] == 10);
        CHECK(y[1] == 20);
        CHECK(y[2] == 30);
        CHECK(y[3] == 40);
        CHECK(y[4] == 50);
    }
}

// TEST_CASE("dataset_basics", "[histogram]") {
//     vector<double> x = {1, 2, 3, 4, 5};
//     vector<double> y = {10, 20, 30, 40, 50};
//     Dataset data(x, y, "i", "j");

//     SECTION("get") {
//         vector<double> i = data.get("i");
//         vector<double> j = data.get("j");
//         CHECK(i == x);
//         CHECK(j == y);
//     }

//     SECTION("reduce") {
//         data.reduce(2);
//         CHECK(data.size() < x.size());
//     }

//     SECTION("limit") {
//         data.limit(Limit(2, 3));
//         vector<double> i = data.get("i");
//         vector<double> j = data.get("j");
//         CHECK(i == vector<double>{2, 3});
//         CHECK(j == vector<double>{20, 30});
//     }
// }

TEST_CASE("dataset_reduce", "[histogram],[files],[manual]") {
    Protein protein("data/2epe.pdb");
    auto h = protein.get_histogram();

    plots::PlotIntensity plot(h);
    SAXSDataset data = h.calc_debye_scattering_intensity();
    data.reduce(20);
    plot.plot_intensity(data);
    plot.save("reduce_test.pdf");
}