#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <settings/SettingsHelper.h>
#include <settings/HistogramSettings.h>
#include <settings/Flags.h>
#include <constants/ConstantsAxes.h>
#include <utility/Utility.h>

using namespace ausaxs;

TEST_CASE("Setting<T>::on_change_and_assignment") {
	SECTION("on_change is called and can modify assigned value") {
		settings::detail::Setting<int> s{0, nullptr};
		bool called = false;

		// on_change will clamp the value to [0, 10] and mark called
		s.on_change = [&](int& v) {
			called = true;
            v = std::clamp(v, 0, 10);
		};

		s = 20; // should be clamped to 10
		CHECK(called == true);
		CHECK(static_cast<int>(s) == 10);

		// assigning negative value should clamp to 0
		called = false;
		s = -5;
		CHECK(called == true);
		CHECK(static_cast<int>(s) == 0);
	}

	SECTION("operator= returns reference to stored value") {
		settings::detail::Setting<int> s{5, nullptr};
		int& r = (s = 42);
		r = 7;
		CHECK(static_cast<int>(s) == 7);
	}
}

TEST_CASE("Setting<T>::conversion_operator") {
	settings::detail::Setting<double> sd{3.14, nullptr};
	double x = sd; // conversion operator
    CHECK_THAT(x, Catch::Matchers::WithinAbs(3.14, 1e-9));
}

TEST_CASE("HistogramSettings::axes::bin_width updates flags") {
	SECTION("setting a custom bin width toggles custom_bin_width and updates inv_bin_width") {
		const double new_width = constants::axes::d_axis.width() * 2.0; // different from default
		settings::axes::bin_width = new_width;
		CHECK(settings::flags::custom_bin_width == true);
        CHECK_THAT(settings::flags::inv_bin_width, Catch::Matchers::WithinAbs(1./new_width, 1e-9));
	}

	SECTION("setting the default width clears custom_bin_width") {
		settings::axes::bin_width = constants::axes::d_axis.width();
		CHECK(settings::flags::custom_bin_width == false);
        CHECK_THAT(settings::flags::inv_bin_width, Catch::Matchers::WithinAbs(1./constants::axes::d_axis.width(), 1e-9));
	}
}