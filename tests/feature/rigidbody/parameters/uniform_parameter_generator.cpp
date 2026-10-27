#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <data/symmetry/BodySymmetryFacade.h>
#include <data/symmetry/CyclicSymmetry.h>
#include <rigidbody/BodySplitter.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <rigidbody/parameters/UniformParameterGenerator.h>
#include <rigidbody/transform/TransformFactory.h>
#include <rigidbody/transform/TransformStrategy.h>  // IWYU pragma: keep
#include <settings/All.h>

#include <rigidbody/parameters/decay/NoDecay.h>

#include <algorithm>
#include <span>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::data;

TEST_CASE("UniformParameterGenerator::next") {
    settings::general::verbose = false;

    int iterations = 100;
    double length_start = GENERATE(1, 2, 3);
    double rad_start = GENERATE(1, 2, 3);
    Rigidbody rb(Molecule{std::vector<Body>{Body(std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)})}});

    SECTION("RotationsOnly") {
        rigidbody::parameter::UniformParameterGenerator ro(&rb, iterations, {.rotation = rad_start});

        for (int i = 0; i < iterations; i++) {
            auto p = ro.next(0);
            CHECK_FALSE(p.translation.has_value());
            // the amplitude bounds the whole rotation, not each of its coordinates
            REQUIRE(p.rotation.value().magnitude() <= rad_start);
        }
    }

    SECTION("TranslationsOnly") {
        rigidbody::parameter::UniformParameterGenerator to(&rb, iterations, {.translation = length_start});

        for (int i = 0; i < iterations; i++) {
            auto p = to.next(0);
            REQUIRE(p.translation.value().magnitude() <= length_start);
            CHECK_FALSE(p.rotation.has_value());
        }
    }

    SECTION("translations & rotations") {
        rigidbody::parameter::UniformParameterGenerator ap(&rb, iterations, {.translation = length_start, .rotation = rad_start});

        for (int i = 0; i < iterations; i++) {
            auto p = ap.next(0);
            REQUIRE(p.translation.value().magnitude() <= length_start);
            REQUIRE(p.rotation.value().magnitude() <= rad_start);
        }
    }
}

TEST_CASE("UniformParameterGenerator::next steps are isotropic") {
    settings::general::verbose = false;
    Rigidbody rb(Molecule{std::vector<Body>{Body(std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)})}});

    // no decay, so that every step is drawn against the full amplitude
    rigidbody::parameter::UniformParameterGenerator gen(
        &rb, std::make_unique<rigidbody::parameter::decay::NoDecay>(), rigidbody::parameter::ParameterAmplitudes{.translation = 1}
    );

    int n = 20000;
    std::vector<int> octants(8, 0);
    Vector3<double> sum = {0, 0, 0};
    double longest = 0;
    for (int i = 0; i < n; ++i) {
        auto t = gen.next(0).translation.value();
        octants[static_cast<int>(t.x() < 0) + 2*static_cast<int>(t.y() < 0) + 4*static_cast<int>(t.z() < 0)]++;
        sum += t/t.magnitude();
        longest = std::max(longest, t.magnitude());
    }

    // every octant is reached about equally often; a cube-shaped draw would pass this too, so the corner test below is the real one
    for (int count : octants) {CHECK(std::abs(count - n/8) < n/20);}

    // the directions cancel out, as they must if no direction is favoured
    CHECK(sum.magnitude()/n < 0.05);

    // the amplitude is a genuine bound: drawing each coordinate independently would reach sqrt(3) beyond it
    CHECK(longest <= 1);
    CHECK(0.99 < longest); // ...and the bound is attained, so it is not merely a loose one
}

TEST_CASE("UniformParameterGenerator::next restricted to a subset of axes") {
    settings::general::verbose = false;
    Rigidbody rb(Molecule{std::vector<Body>{Body(std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)})}});

    SECTION("translations in a plane stay in it, and are isotropic within it") {
        rigidbody::parameter::UniformParameterGenerator gen(
            &rb, std::make_unique<rigidbody::parameter::decay::NoDecay>(),
            rigidbody::parameter::ParameterAmplitudes{.translation = 1, .translation_axes = {.x = true, .y = true, .z = false}}
        );

        int n = 20000;
        std::vector<int> quadrants(4, 0);
        double longest = 0;
        for (int i = 0; i < n; ++i) {
            auto t = gen.next(0).translation.value();
            REQUIRE(t.z() == 0);
            quadrants[static_cast<int>(t.x() < 0) + 2*static_cast<int>(t.y() < 0)]++;
            longest = std::max(longest, t.magnitude());
        }
        for (int count : quadrants) {CHECK(std::abs(count - n/4) < n/20);}
        CHECK(longest <= 1);
        CHECK(0.99 < longest);
    }

    SECTION("a single axis draws both directions along it") {
        rigidbody::parameter::UniformParameterGenerator gen(
            &rb, std::make_unique<rigidbody::parameter::decay::NoDecay>(),
            rigidbody::parameter::ParameterAmplitudes{.rotation = 1, .rotation_axes = {.x = false, .y = false, .z = true}}
        );

        int n = 2000, positive = 0;
        for (int i = 0; i < n; ++i) {
            auto r = gen.next(0).rotation.value();
            REQUIRE(r.x() == 0);
            REQUIRE(r.y() == 0);
            REQUIRE(std::abs(r.z()) <= 1);
            if (0 < r.z()) {++positive;}
        }
        CHECK(std::abs(positive - n/2) < n/10);
    }
}

namespace {
    // the three ways a body move is carried out: about a constraint pivot for the body alone, about a constraint pivot for the
    // whole branch connected through it, or about the body's own centre of mass when it is moved without a constraint
    enum class MovePath {single, rigid, free};

    struct PlanarFixture {
        PlanarFixture(MovePath path) : path(path), rb(BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99})) {
            rb.transformer = factory::create_transform_strategy(&rb, path == MovePath::single 
                ? settings::rigidbody::TransformationStrategyChoice::SingleTransform 
                : settings::rigidbody::TransformationStrategyChoice::RigidTransform
            );
            rb.constraints->generate_constraints(settings::rigidbody::ConstraintGenerationStrategyChoice::Backbone);
            rb.molecule.generate_new_hydration();
        }

        void use(const rigidbody::parameter::ParameterAmplitudes& amplitudes) {
            rb.parameter_generator = std::make_shared<rigidbody::parameter::UniformParameterGenerator>(
                &rb, std::make_unique<rigidbody::parameter::decay::NoDecay>(), amplitudes
            );
        }

        void step(int ibody) {
            auto param = rb.parameter_generator->next(ibody);
            if (path == MovePath::free) {rb.transformer->apply(param, ibody);}
            else {rb.transformer->apply(param, rb.constraints->get_body_constraints(ibody).at(0), ibody);}
        }

        void steps(int n) {for (int i = 0; i < n; ++i) {step(i % rb.molecule.size_body());}}

        std::vector<double> cm_heights() {
            std::vector<double> out;
            out.reserve(rb.molecule.size_body());
            for (int i = 0; i < rb.molecule.size_body(); ++i) {out.push_back(rb.molecule.get_body(i).get_cm().z());}
            return out;
        }

        std::vector<double> atom_heights() {
            std::vector<double> out;
            for (int i = 0; i < rb.molecule.size_body(); ++i) {
                for (const auto& a : rb.molecule.get_body(i).get_atoms()) {out.push_back(a.coordinates().z());}
            }
            return out;
        }

        static double largest_change(const std::vector<double>& before, const std::vector<double>& after) {
            REQUIRE(before.size() == after.size());
            double worst = 0;
            for (std::size_t i = 0; i < before.size(); ++i) {worst = std::max(worst, std::abs(after[i] - before[i]));}
            return worst;
        }

        MovePath path;
        Rigidbody rb;
    };
}

TEST_CASE("UniformParameterGenerator::next in-plane translations keep every unconstrained body centred in its plane, however it rotates") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::grid::min_bins = 250;

    // only a body moved on its own is rotated about its centre of mass; a constrained move rotates about the constraint pivot and so
    // displaces the centre of mass, which is inherent to it
    PlanarFixture f(MovePath::free);
    f.use({.translation = 3, .rotation = 1, .translation_axes = {.x = true, .y = true, .z = false}});

    auto initial_cm = f.cm_heights();
    auto initial_atoms = f.atom_heights();
    f.steps(30);

    REQUIRE_THAT(PlanarFixture::largest_change(initial_cm, f.cm_heights()), Catch::Matchers::WithinAbs(0, 1e-9));
    CHECK(1 < PlanarFixture::largest_change(initial_atoms, f.atom_heights())); // ...while the bodies did tilt out of the plane
}

TEST_CASE("UniformParameterGenerator::next rotations about z leave every atom at its height, also for a tilted body") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::grid::min_bins = 250;

    PlanarFixture f(GENERATE(MovePath::single, MovePath::rigid, MovePath::free));

    // tilt the bodies first, so that their own z-axes no longer coincide with the z-axis of the frame they were loaded in
    f.use({.rotation = 1});
    f.steps(6);

    f.use({.translation = 3, .rotation = 1, .translation_axes = {.x = true, .y = true, .z = false}, .rotation_axes = {.x = false, .y = false, .z = true}});
    auto initial = f.atom_heights();
    f.steps(30);
    REQUIRE_THAT(PlanarFixture::largest_change(initial, f.atom_heights()), Catch::Matchers::WithinAbs(0, 1e-9));
}

TEST_CASE("UniformParameterGenerator::next symmetry components are independent") {
    settings::general::verbose = false;

    int iterations = 100;
    Molecule m{std::vector<Body>{Body(std::vector{
        AtomFF({0, 0, 0}, form_factor::form_factor_t::C), AtomFF({1, 1, 1}, form_factor::form_factor_t::C)
    })}};
    m.get_body(0).symmetry().add(symmetry::type::c2);
    Rigidbody rb(std::move(m));

    auto all_zero = [](std::span<double> s) {return std::ranges::all_of(s, [](double v) {return v == 0;});};
    auto any_nonzero = [](std::span<double> s) {return std::ranges::any_of(s, [](double v) {return v != 0;});};

    SECTION("only the symmetry translation moves when only its amplitude is set") {
        rigidbody::parameter::UniformParameterGenerator gen(&rb, iterations, {.symmetry_translation = 5});

        bool saw_translation = false;
        for (int i = 0; i < iterations; i++) {
            auto p = gen.next(0);
            CHECK_FALSE(p.translation.has_value());
            CHECK_FALSE(p.rotation.has_value());
            REQUIRE(p.symmetry_pars.has_value());
            REQUIRE(p.symmetry_pars.value().size() == 1);
            auto& delta = *p.symmetry_pars.value()[0];
            CHECK(all_zero(delta.span_rotation()));
            saw_translation |= any_nonzero(delta.span_translation());
        }
        CHECK(saw_translation);
    }

    SECTION("only the symmetry rotation moves when only its amplitude is set") {
        rigidbody::parameter::UniformParameterGenerator gen(&rb, iterations, {.symmetry_rotation = 0.5});

        bool saw_rotation = false;
        for (int i = 0; i < iterations; i++) {
            auto p = gen.next(0);
            CHECK_FALSE(p.translation.has_value());
            CHECK_FALSE(p.rotation.has_value());
            REQUIRE(p.symmetry_pars.has_value());
            REQUIRE(p.symmetry_pars.value().size() == 1);
            auto& delta = *p.symmetry_pars.value()[0];
            CHECK(all_zero(delta.span_translation()));
            saw_rotation |= any_nonzero(delta.span_rotation());
        }
        CHECK(saw_rotation);
    }

    SECTION("no symmetry parameters are generated when neither amplitude is set") {
        rigidbody::parameter::UniformParameterGenerator gen(&rb, iterations, {.translation = 1, .rotation = 0.1});
        for (int i = 0; i < iterations; i++) {
            auto p = gen.next(0);
            CHECK(p.translation.has_value());
            CHECK(p.rotation.has_value());
            CHECK_FALSE(p.symmetry_pars.has_value());
        }
    }
}

TEST_CASE("UniformParameterGenerator: a cyclic axis keeps its length as deltas accumulate") {
    settings::general::verbose = false;

    Molecule m{std::vector<Body>{Body(std::vector{
        AtomFF({0, 0, 0}, form_factor::form_factor_t::C), AtomFF({1, 1, 1}, form_factor::form_factor_t::C)
    })}};
    m.get_body(0).symmetry().add(symmetry::type::c2);
    Rigidbody rb(std::move(m));

    // the axis is a direction, so only its orientation is meaningful. Before the amplitude was an angle, deltas were added
    // to it component-wise and its length random-walked upward, quietly shrinking the rotation a given amplitude produced.
    rb.parameter_generator = std::make_shared<rigidbody::parameter::UniformParameterGenerator>(
        &rb, 100, rigidbody::parameter::ParameterAmplitudes{.symmetry_rotation = 0.2}
    );

    for (int i = 0; i < 50; ++i) {
        rb.transformer->apply(rb.parameter_generator->next(0), 0);
        auto* sym = static_cast<symmetry::CyclicSymmetry*>(rb.molecule.get_body(0).symmetry().get(0));
        REQUIRE_THAT(sym->_repeat_relation.axis.magnitude(), Catch::Matchers::WithinAbs(1, 1e-9));
    }
}

TEST_CASE("UniformParameterGenerator::next respects a planar symmetry's reduced parameter count") {
    settings::general::verbose = false;

    Molecule m{std::vector<Body>{Body(std::vector{
        AtomFF({0, 0, 0}, form_factor::form_factor_t::C), AtomFF({1, 1, 1}, form_factor::form_factor_t::C)
    })}};
    m.get_body(0).symmetry().add(symmetry::type::dp3); // a planar dihedral exposes only 2 of its 3 translation components
    Rigidbody rb(std::move(m));

    rigidbody::parameter::UniformParameterGenerator gen(
        &rb, std::make_unique<rigidbody::parameter::decay::NoDecay>(), rigidbody::parameter::ParameterAmplitudes{.symmetry_translation = 4}
    );

    double longest = 0;
    for (int i = 0; i < 2000; ++i) {
        auto p = gen.next(0);
        auto t = p.symmetry_pars.value()[0]->span_translation();
        REQUIRE(t.size() == 2);
        double magnitude = std::sqrt(t[0]*t[0] + t[1]*t[1]);

        // the direction is drawn within the plane the symmetry actually exposes, so the amplitude bounds the step there too.
        // projecting a three-dimensional direction onto the plane would instead shorten it by a factor of sin(polar angle).
        REQUIRE(magnitude <= 4);
        longest = std::max(longest, magnitude);
    }
    CHECK(3.99 < longest);
}

TEST_CASE("UniformParameterGenerator: a rejected step does not corrupt the next delta") {
    settings::general::verbose = false;

    Molecule m{std::vector<Body>{Body(std::vector{
        AtomFF({0, 0, 0}, form_factor::form_factor_t::C), AtomFF({1, 1, 1}, form_factor::form_factor_t::C)
    })}};
    m.get_body(0).symmetry().add(symmetry::type::c2);
    Rigidbody rb(std::move(m));

    rb.parameter_generator = std::make_shared<rigidbody::parameter::UniformParameterGenerator>(
        &rb, std::make_unique<rigidbody::parameter::decay::NoDecay>(), rigidbody::parameter::ParameterAmplitudes{.symmetry_rotation = 0.2}
    );

    // a symmetry-only apply consumes the body backup, so undo restores the parameters but leaves the body holding the
    // rejected axis. The generator must read the parameters rather than the body, or the delta is computed against the
    // rejected axis and added to the restored one, which neither preserves the unit length nor the requested angle.
    for (int i = 0; i < 50; ++i) {
        rb.transformer->apply(rb.parameter_generator->next(0), 0);
        rb.transformer->undo();

        rb.transformer->apply(rb.parameter_generator->next(0), 0);
        const auto& params = rb.conformation->absolute_parameters.parameters[0].symmetry_pars;
        auto* sym = static_cast<symmetry::CyclicSymmetry*>(params[0].get());
        REQUIRE_THAT(sym->_repeat_relation.axis.magnitude(), Catch::Matchers::WithinAbs(1, 1e-9));
    }
}
