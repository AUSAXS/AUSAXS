#include <catch2/catch_test_macros.hpp>

#include <plots/PlotOptions.h>

using namespace ausaxs;

TEST_CASE("PlotOptions::set") {
    plots::PlotOptions options;
    options.set({{"xlabel", "label for x"}, {"ylabel", std::string("label for y")}, {"line", true}, {"errors", true}, {"color", style::color::blue}, {"alpha", 0.3}});
    CHECK(options.xlabel == "label for x");
    CHECK(options.ylabel == "label for y");
    CHECK(options.draw_line == true);
    CHECK(options.draw_errors == true);
    CHECK(options.color == style::color::blue);
    CHECK(options.alpha == 0.3);
}
