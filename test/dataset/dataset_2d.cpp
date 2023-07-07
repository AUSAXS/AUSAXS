TEST_CASE("Dataset2D_scaling_methods") {
    std::vector<double> x = {1, 2, 3, 4, 5};
    std::vector<double> y = {10, 20, 30, 40, 50};
    std::vector<double> yerr = {1, 2, 3, 4, 5};
    vector<double> xerr = {0.1, 0.2, 0.3, 0.4, 0.5};
    Dataset2D data(x, y, xerr, yerr);

    SECTION("scale y") {
        data.scale_y(2);
        CHECK(data.y() == std::vector<double>({20, 40, 60, 80, 100}));

        data.scale_y(0.2);
        CHECK(data.y() == std::vector<double>({4, 8, 12, 16, 20}));

        data.scale_y(-0.5);
        CHECK(data.y() == std::vector<double>({-2, -4, -6, -8, -10}));
    }

    SECTION("scale errors") {
        data.scale_errors(2);
        CHECK(data.yerr() == std::vector<double>({2, 4, 6, 8, 10}));
        CHECK(data.xerr() == std::vector<double>({0.2, 0.4, 0.6, 0.8, 1.0}));

        data.scale_errors(0.5);
        CHECK(data.yerr() == std::vector<double>({1, 2, 3, 4, 5}));
        CHECK(data.xerr() == std::vector<double>({0.1, 0.2, 0.3, 0.4, 0.5}));

        data.scale_errors(-0.5);
        CHECK(data.yerr() == std::vector<double>({-0.5, -1, -1.5, -2, -2.5}));
        CHECK(data.xerr() == std::vector<double>({-0.05, -0.1, -0.15, -0.2, -0.25}));
    }

    SECTION("normalize") {
        data.normalize(1);
        CHECK(data.y() == std::vector<double>({1, 2, 3, 4, 5}));

        data.normalize(-5);
        CHECK(data.y() == std::vector<double>({-5, -10, -15, -20, -25}));

        data.normalize(10);
        CHECK(data.y() == std::vector<double>({10, 20, 30, 40, 50}));
    }
}