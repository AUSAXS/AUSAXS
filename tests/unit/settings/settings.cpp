#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <constants/ConstantsAxes.h>
#include <settings/HistogramSettings.h>
#include <settings/InternalState.h>
#include <settings/SettingsHelper.h>

#include <type_traits>

using namespace ausaxs;

TEST_CASE("Setting<T>::on_change_and_assignment") {
	SECTION("on_change is called and can modify assigned value") {
		settings::detail::Setting<int> s{0};
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
		settings::detail::Setting<int> s{5};
		int& r = (s = 42);
		r = 7;
		CHECK(static_cast<int>(s) == 7);
	}
}

TEST_CASE("Setting<T>::conversion_operator") {
	settings::detail::Setting<double> sd{3.14};
	double x = sd; // conversion operator
    CHECK_THAT(x, Catch::Matchers::WithinAbs(3.14, 1e-9));
}

TEST_CASE("Setting<T>::not_copyable") {
	// copying a setting and assigning it back would silently skip on_change
	STATIC_REQUIRE_FALSE(std::is_copy_constructible_v<settings::detail::Setting<int>>);
	STATIC_REQUIRE_FALSE(std::is_copy_assignable_v<settings::detail::Setting<int>>);
	STATIC_REQUIRE(std::is_assignable_v<settings::detail::Setting<int>&, int>);
}

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