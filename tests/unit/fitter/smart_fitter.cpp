#include <catch2/catch_test_macros.hpp>

#include <constants/ConstantsFitParameters.h>
#include <dataset/SimpleDataset.h>
#include <fitter/SmartFitter.h>
#include <mini/detail/Parameter.h>
#include <utility/Exceptions.h>

using namespace ausaxs;
using namespace ausaxs::fitter;

namespace {
    // exposes the reordered guess so the permutation performed by set_guess can be inspected
    struct GuessProbe : public SmartFitter {
        using SmartFitter::SmartFitter;
        const std::vector<mini::Parameter>& get_guess() const {return guess;}
    };

    std::string name_of(constants::fit::Parameters p) {return constants::fit::to_string(p);}

    GuessProbe make_probe(SmartFitter::EnabledFitParameters enabled) {
        GuessProbe probe(SimpleDataset({1, 2, 3}, {1, 1, 1}, {1, 1, 1}));
        probe.enabled_fit_parameters = enabled;
        return probe;
    }
}

TEST_CASE("SmartFitter::set_guess") {
    auto cw = name_of(constants::fit::Parameters::SCALING_WATER);
    auto cx = name_of(constants::fit::Parameters::SCALING_EXV);
    auto cr = name_of(constants::fit::Parameters::SCALING_RHO);

    SECTION("reorders into canonical order") {
        auto probe = make_probe({.hydration=true, .excluded_volume=true, .solvent_density=true, .atomic_debye_waller=false, .exv_debye_waller=false});
        probe.set_guess({
            mini::Parameter{cr, 3, {3.1, 3.2}},
            mini::Parameter{cw, 1, {1.1, 1.2}},
            mini::Parameter{cx, 2, {2.1, 2.2}}
        });

        const auto& g = probe.get_guess();
        REQUIRE(g.size() == 3);

        // each parameter must keep its own guess *and* its own bounds
        CHECK(g[0].name == cw);
        CHECK(g[0].guess.value() == 1);
        CHECK(g[0].bounds.value().min == 1.1);
        CHECK(g[1].name == cx);
        CHECK(g[1].guess.value() == 2);
        CHECK(g[1].bounds.value().min == 2.1);
        CHECK(g[2].name == cr);
        CHECK(g[2].guess.value() == 3);
        CHECK(g[2].bounds.value().min == 3.1);
    }

    SECTION("handles non-prefix parameter sets") {
        // hydration disabled, so the enabled slots are {1, 2} rather than {0, 1}
        auto probe = make_probe({.hydration=false, .excluded_volume=true, .solvent_density=true, .atomic_debye_waller=false, .exv_debye_waller=false});
        probe.set_guess({mini::Parameter{cx, 2}, mini::Parameter{cr, 3}});

        const auto& g = probe.get_guess();
        REQUIRE(g.size() == 2);
        CHECK(g[0].name == cx);
        CHECK(g[0].guess.value() == 2);
        CHECK(g[1].name == cr);
        CHECK(g[1].guess.value() == 3);
    }

    SECTION("rejects disabled and unknown parameters") {
        auto probe = make_probe({.hydration=true, .excluded_volume=false, .solvent_density=false, .atomic_debye_waller=false, .exv_debye_waller=false});
        CHECK_THROWS_AS(probe.set_guess({mini::Parameter{cx, 1}}), except::invalid_argument);
        CHECK_THROWS_AS(probe.set_guess({mini::Parameter{"nonsense", 1}}), except::invalid_argument);
        CHECK_THROWS_AS(probe.set_guess({mini::Parameter{cw, 1}, mini::Parameter{cx, 2}}), except::invalid_argument);
    }
}
