#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <settings/All.h>
#include <settings/SettingRef.h>
#include <settings/SettingsHelper.h>

#include <map>
#include <string>
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

TEST_CASE("Setting<T>::reset") {
	SECTION("restores the constructed value through on_change") {
		settings::detail::Setting<int> s{3};
		int calls = 0;
		s.on_change = [&](int&) {++calls;};
		s = 8;
		s.reset();
		CHECK(static_cast<int>(s) == 3);
		CHECK(calls == 2);
	}

	SECTION("settings::reset restores every global setting") {
		settings::reset();
		std::map<std::string, std::string> defaults;
		for (auto& [name, ref] : settings::io::detail::ISettingRef::get_stored_settings()) {defaults[name] = ref->get();}

		// the writes EM, rigidbody and a previous CLI run leave behind
		settings::exv::exv_method = settings::exv::ExvMethod::None;
		settings::molecule::implicit_hydrogens = false;
		settings::grid::min_exv_radius = 3.7;
		settings::grid::scaling = 2;
		settings::grid::cubic = true;
		settings::axes::qmax = 0.2;
		settings::axes::bin_width = 0.5;
		settings::general::output = "elsewhere/";
		settings::reset();

		for (auto& [name, ref] : settings::io::detail::ISettingRef::get_stored_settings()) {
			INFO(name);
			CHECK(ref->get() == defaults[name]);
		}
		CHECK(settings::internal_state::inv_bin_width == 1/settings::axes::bin_width); // derived state follows through on_change
	}
}
