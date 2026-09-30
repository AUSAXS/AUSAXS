#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <constants/ConstantsAxes.h>
#include <settings/HistogramSettings.h>
#include <settings/InternalState.h>

using namespace ausaxs;

TEST_CASE("HistogramSettings::axes::bin_width updates inv_bin_width") {
	SECTION("setting a custom bin width updates inv_bin_width") {
		const double new_width = constants::axes::d_axis.width() * 2.0; // different from default
		settings::axes::bin_width = new_width;
        CHECK_THAT(settings::internal_state::inv_bin_width, Catch::Matchers::WithinAbs(1./new_width, 1e-9));
	}

	SECTION("setting the default width restores inv_bin_width") {
		settings::axes::bin_width = constants::axes::d_axis.width();
        CHECK_THAT(settings::internal_state::inv_bin_width, Catch::Matchers::WithinAbs(1./constants::axes::d_axis.width(), 1e-9));
	}
}
