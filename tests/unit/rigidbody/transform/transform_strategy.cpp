#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <data/detail/SimpleBody.h>
#include <data/symmetry/CyclicSymmetry.h>
#include <data/symmetry/PointSymmetry.h>
#include <math/MatrixUtils.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/constraints/DistanceConstraintBond.h>
#include <rigidbody/constraints/IDistanceConstraint.h>
#include <rigidbody/detail/MoleculeTransformParametersAbsolute.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <rigidbody/parameters/BodyTransformParametersAbsolute.h>
#include <rigidbody/parameters/BodyTransformParametersRelative.h>
#include <rigidbody/transform/RigidTransform.h>
#include <rigidbody/transform/SingleTransform.h>
#include <rigidbody/transform/TransformGroup.h>
#include <rigidbody/transform/TransformStrategy.h>  // IWYU pragma: keep
#include <settings/All.h>

#include <support/rb_metadata.h>

#include <algorithm>
#include <numbers>
#include <unordered_set>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;

TEST_CASE("TransformStrategy::apply unconstrained body") {
    settings::molecule::implicit_hydrogens = false;
    settings::grid::min_bins = 25;

    SECTION("translate unconstrained body") {
        AtomFF a1({1, 2, 3}, form_factor::form_factor_t::C);
        Body b1(std::vector<AtomFF>{a1});
        
        Rigidbody rigidbody(Molecule{std::vector<Body>{b1}});
        auto& transformer = rigidbody.transformer;
        
        auto cm_before = rigidbody.molecule.get_body(0).get_cm();
        
        transformer->apply({{1, 1, 1}, {0, 0, 0}}, 0u);
        
        auto cm_after = rigidbody.molecule.get_body(0).get_cm();
        REQUIRE_THAT(cm_after.x(), Catch::Matchers::WithinAbs(cm_before.x() + 1, 1e-10));
        REQUIRE_THAT(cm_after.y(), Catch::Matchers::WithinAbs(cm_before.y() + 1, 1e-10));
        REQUIRE_THAT(cm_after.z(), Catch::Matchers::WithinAbs(cm_before.z() + 1, 1e-10));
    }

    SECTION("rotate unconstrained body") {
        AtomFF a1({1, 0, 0}, form_factor::form_factor_t::C);
        Body b1(std::vector<AtomFF>{a1});
        
        Rigidbody rigidbody(Molecule{std::vector<Body>{b1}});
        auto& transformer = rigidbody.transformer;
        
        auto cm_before = rigidbody.molecule.get_body(0).get_cm();
        
        transformer->apply({{0, 0, 0}, {0, 0, std::numbers::pi/2}}, 0u);
        
        auto cm = rigidbody.molecule.get_body(0).get_cm();
        // Since the body is centered in initial_conformation at (0,0,0),
        // rotating it around its own center doesn't change its position
        // The translation component of absolute_parameters keeps it at the same place
        REQUIRE_THAT(cm.x(), Catch::Matchers::WithinAbs(cm_before.x(), 1e-10));
        REQUIRE_THAT(cm.y(), Catch::Matchers::WithinAbs(cm_before.y(), 1e-10));
    }

    SECTION("parameters updated for unconstrained body") {
        AtomFF a1({5, 5, 5}, form_factor::form_factor_t::C);
        Body b1(std::vector<AtomFF>{a1});
        
        Rigidbody rigidbody(Molecule{std::vector<Body>{b1}});
        auto& transformer = rigidbody.transformer;
        
        auto params_before = rigidbody.conformation->absolute_parameters.parameters[0];
        
        transformer->apply({{2, 3, 4}, {0.1, 0.2, 0.3}}, 0u);
        
        auto params_after = rigidbody.conformation->absolute_parameters.parameters[0];
        REQUIRE(params_after.translation != params_before.translation);
        REQUIRE(params_after.rotation != params_before.rotation);
    }
}

TEST_CASE("TransformStrategy::rotate_and_translate") {
    settings::molecule::implicit_hydrogens = false;
    settings::grid::min_bins = 10;

    SECTION("rotation then translation on single body") {
        AtomFF a1({1, 0, 0}, form_factor::form_factor_t::C);
        Body b1(std::vector<AtomFF>{a1});
        
        Rigidbody rigidbody(Molecule{std::vector<Body>{b1}});
        
        auto& body = rigidbody.molecule.get_body(0);
        auto cm_before = body.get_cm();
        
        // Apply via the public interface that uses rotate_and_translate
        auto& transformer = rigidbody.transformer;
        transformer->apply({{1, 1, 0}, {0, 0, std::numbers::pi/2}}, 0u);
        
        auto cm = body.get_cm();
        // Body is centered in initial_conformation, so rotation doesn't change position
        // Translation adds (1,1,0) to the current position
        REQUIRE_THAT(cm.x(), Catch::Matchers::WithinAbs(cm_before.x() + 1.0, 1e-10));
        REQUIRE_THAT(cm.y(), Catch::Matchers::WithinAbs(cm_before.y() + 1.0, 1e-10));
    }
}

TEST_CASE("TransformStrategy::undo") {
    settings::molecule::implicit_hydrogens = false;
    settings::grid::min_bins = 50;

    SECTION("undo after single transformation") {
        AtomFF a1({3, 4, 5}, form_factor::form_factor_t::C);
        Body b1(std::vector<AtomFF>{a1});
        
        Rigidbody rigidbody(Molecule{std::vector<Body>{b1}});
        auto& transformer = rigidbody.transformer;
        
        auto cm_original = rigidbody.molecule.get_body(0).get_cm();
        auto params_original = rigidbody.conformation->absolute_parameters.parameters[0];
        
        transformer->apply({{5, 6, 7}, {1, 2, 3}}, 0u);
        
        // Verify change
        auto cm_after = rigidbody.molecule.get_body(0).get_cm();
        REQUIRE(cm_after != cm_original);
        
        // Undo
        transformer->undo();
        
        // Verify restoration
        auto cm_restored = rigidbody.molecule.get_body(0).get_cm();
        auto params_restored = rigidbody.conformation->absolute_parameters.parameters[0];
        
        REQUIRE_THAT(cm_restored.x(), Catch::Matchers::WithinAbs(cm_original.x(), 1e-10));
        REQUIRE_THAT(cm_restored.y(), Catch::Matchers::WithinAbs(cm_original.y(), 1e-10));
        REQUIRE_THAT(cm_restored.z(), Catch::Matchers::WithinAbs(cm_original.z(), 1e-10));
        
        REQUIRE_THAT(params_restored.translation.x(), Catch::Matchers::WithinAbs(params_original.translation.x(), 1e-10));
        REQUIRE_THAT(params_restored.rotation.x(), Catch::Matchers::WithinAbs(params_original.rotation.x(), 1e-10));
    }

    SECTION("multiple undos only restore last transformation") {
        AtomFF a1({0, 0, 0}, form_factor::form_factor_t::C);
        Body b1(std::vector<AtomFF>{a1});
        
        Rigidbody rigidbody(Molecule{std::vector<Body>{b1}});
        auto& transformer = rigidbody.transformer;
        
        transformer->apply({{1, 0, 0}, {0, 0, 0}}, 0u);
        auto cm_after_first = rigidbody.molecule.get_body(0).get_cm();
        
        transformer->apply({{1, 0, 0}, {0, 0, 0}}, 0u);
        // auto cm_after_second = rigidbody.molecule.get_body(0).get_cm();
        
        // Undo should restore to state after first transform
        transformer->undo();
        auto cm_after_undo = rigidbody.molecule.get_body(0).get_cm();
        
        REQUIRE_THAT(cm_after_undo.x(), Catch::Matchers::WithinAbs(cm_after_first.x(), 1e-10));
        REQUIRE_THAT(cm_after_undo.y(), Catch::Matchers::WithinAbs(cm_after_first.y(), 1e-10));
    }
}

TEST_CASE("TransformStrategy::reconstructed body matches current state after multiple transformations") {
    settings::molecule::implicit_hydrogens = false;
    settings::grid::min_bins = 10;

    SECTION("reconstructed body matches current state after multiple transformations") {
        AtomFF a1({2, 3, 4}, form_factor::form_factor_t::C);
        AtomFF a2({5, 6, 7}, form_factor::form_factor_t::C);
        Body b1(std::vector<AtomFF>{a1, a2});
        
        Rigidbody rigidbody(Molecule{std::vector<Body>{b1}});
        auto& transformer = rigidbody.transformer;
        
        // Apply multiple transformations
        transformer->apply({{1, 2, 3}, {0.5, 0.3, 0.1}}, 0u);
        transformer->apply({{-1, 1, -1}, {0.2, -0.1, 0.4}}, 0u);
        
        // Reconstruct from initial conformation
        auto& current_body = rigidbody.molecule.get_body(0);
        auto& params = rigidbody.conformation->absolute_parameters.parameters[0];
        auto& initial = rigidbody.conformation->initial_conformation[0];
        
        Body reconstructed = initial;
        reconstructed.rotate(matrix::rotation_matrix(params.rotation));
        reconstructed.translate(params.translation);
        
        // Compare all atoms
        REQUIRE(reconstructed.size_atom() == current_body.size_atom());
        for (int i = 0; i < current_body.size_atom(); ++i) {
            auto current_pos = current_body.get_atom(i).coordinates();
            auto reconstructed_pos = reconstructed.get_atom(i).coordinates();
            
            REQUIRE_THAT(reconstructed_pos.x(), Catch::Matchers::WithinAbs(current_pos.x(), 1e-3));
            REQUIRE_THAT(reconstructed_pos.y(), Catch::Matchers::WithinAbs(current_pos.y(), 1e-3));
            REQUIRE_THAT(reconstructed_pos.z(), Catch::Matchers::WithinAbs(current_pos.z(), 1e-3));
        }
    }
}

TEST_CASE("TransformStrategy::parameter accumulation") {
    settings::molecule::implicit_hydrogens = false;
    settings::grid::min_bins = 10;

    SECTION("translations accumulate linearly") {
        AtomFF a1({0, 0, 0}, form_factor::form_factor_t::C);
        Body b1(std::vector<AtomFF>{a1});
        
        Rigidbody rigidbody(Molecule{std::vector<Body>{b1}});
        auto& transformer = rigidbody.transformer;
        
        transformer->apply({{1, 0, 0}, {0, 0, 0}}, 0u);
        transformer->apply({{0, 2, 0}, {0, 0, 0}}, 0u);
        transformer->apply({{0, 0, 3}, {0, 0, 0}}, 0u);
        
        auto& params = rigidbody.conformation->absolute_parameters.parameters[0];
        REQUIRE_THAT(params.translation.x(), Catch::Matchers::WithinAbs(1.0, 1e-10));
        REQUIRE_THAT(params.translation.y(), Catch::Matchers::WithinAbs(2.0, 1e-10));
        REQUIRE_THAT(params.translation.z(), Catch::Matchers::WithinAbs(3.0, 1e-10));
    }

    SECTION("rotations accumulate correctly") {
        AtomFF a1({1, 0, 0}, form_factor::form_factor_t::C);
        Body b1(std::vector<AtomFF>{a1});
        
        Rigidbody rigidbody(Molecule{std::vector<Body>{b1}});
        auto& transformer = rigidbody.transformer;
        
        auto cm_before = rigidbody.molecule.get_body(0).get_cm();
        
        // Apply four 45-degree rotations = 180 degrees total
        for (int i = 0; i < 4; ++i) {
            transformer->apply({{0, 0, 0}, {0, 0, std::numbers::pi/4}}, 0u);
        }
        
        auto& params = rigidbody.conformation->absolute_parameters.parameters[0];
        REQUIRE_THAT(params.rotation.z(), Catch::Matchers::WithinAbs(std::numbers::pi, 1e-6));
        
        // Position should remain the same since body is centered in initial_conformation
        // and we're only applying rotations (no translation)
        auto cm = rigidbody.molecule.get_body(0).get_cm();
        REQUIRE_THAT(cm.x(), Catch::Matchers::WithinAbs(cm_before.x(), 1e-10));
        REQUIRE_THAT(cm.y(), Catch::Matchers::WithinAbs(cm_before.y(), 1e-10));
    }
}

// A symmetry is anchored to the body it belongs to, so a real-space move must carry its copies along rigidly. The assembly is then invariant
// under everything the optimiser can do in real space, leaving the symmetry parameters as the only ones that can change its shape.
TEST_CASE("TransformStrategy::apply keeps a symmetric assembly rigid") {
    settings::molecule::implicit_hydrogens = false;
    settings::grid::min_bins = 25;

    auto make_body = [] (std::unique_ptr<symmetry::ISymmetry> sym) {
        Body b(std::vector<AtomFF>{
            AtomFF({ 1.0,  0.0,  0.0}, form_factor::form_factor_t::C),
            AtomFF({-1.0,  0.5,  0.2}, form_factor::form_factor_t::C),
            AtomFF({ 0.3, -1.4,  0.9}, form_factor::form_factor_t::C),
            AtomFF({ 2.1,  0.7, -1.1}, form_factor::form_factor_t::C)
        });
        b.symmetry().add(std::move(sym));
        return b;
    };

    // every pairwise distance in the full assembly; invariant under any rigid motion of it, and not under a deformation
    auto assembly_distances = [] (const Rigidbody& rigidbody) {
        auto structure = rigidbody.molecule.get_body(0).symmetry().explicit_structure();
        std::vector<double> distances;
        for (std::size_t i = 0; i < structure.atoms.size(); ++i) {
            for (std::size_t j = i+1; j < structure.atoms.size(); ++j) {
                distances.push_back(structure.atoms[i].coordinates().distance(structure.atoms[j].coordinates()));
            }
        }
        std::ranges::sort(distances);
        return distances;
    };

    // a real-space-only move: the symmetry parameters are left untouched, so only the base body is driven
    auto real_space_move = [] (const Vector3<double>& translation, const Vector3<double>& rotation) {
        parameter::BodyTransformParametersRelative par;
        par.translation = translation;
        par.rotation = rotation;
        return par;
    };

    auto check_invariant = [&] (Rigidbody& rigidbody) {
        auto before = assembly_distances(rigidbody);
        rigidbody.transformer->apply(real_space_move({3, -2, 5}, {0.4, -0.9, 0.3}), 0u);
        auto after = assembly_distances(rigidbody);

        // a second move, to catch a pose that is accumulated rather than applied absolutely
        rigidbody.transformer->apply(real_space_move({-1, 4, 2}, {-0.7, 0.2, 1.1}), 0u);
        auto after2 = assembly_distances(rigidbody);

        REQUIRE(before.size() == after.size());
        REQUIRE(before.size() == after2.size());
        for (std::size_t i = 0; i < before.size(); ++i) {
            REQUIRE_THAT(after[i],  Catch::Matchers::WithinAbs(before[i], 1e-9));
            REQUIRE_THAT(after2[i], Catch::Matchers::WithinAbs(before[i], 1e-9));
        }
    };

    SECTION("point symmetry") {
        Rigidbody rigidbody(Molecule{std::vector<Body>{
            make_body(std::make_unique<symmetry::PointSymmetry>(Vector3<double>{8, 0, 0}, Vector3<double>{0.3, 0.7, 0.2}))
        }});
        check_invariant(rigidbody);
    }

    SECTION("cyclic symmetry") {
        Rigidbody rigidbody(Molecule{std::vector<Body>{
            make_body(std::make_unique<symmetry::CyclicSymmetry>(
                Vector3<double>{8, 0, 0}, Vector3<double>{0, 0, 0}, Vector3<double>{0, 0, 1}, 2*std::numbers::pi/3, 2
            ))
        }});
        check_invariant(rigidbody);
    }
}

TEST_CASE("TransformStrategy::apply zeroed delta parameters do not cause drift") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;

    AtomFF a1({0, 0, 0}, form_factor::form_factor_t::C);
    AtomFF a2({5, 0, 0}, form_factor::form_factor_t::C);
    AtomFF a3({0, 5, 0}, form_factor::form_factor_t::C);
    Body b1(std::vector<AtomFF>{a1});
    Body b2(std::vector<AtomFF>{a2});
    Body b3(std::vector<AtomFF>{a3});
    Rigidbody rigidbody(Molecule{std::vector<Body>{b1, b2, b3}});

    auto& transformer = rigidbody.transformer;

    SECTION("repeated zero transforms do not drift rotation") {
        auto original_rotation = rigidbody.conformation->absolute_parameters.parameters[0].rotation;
        auto original_translation = rigidbody.conformation->absolute_parameters.parameters[0].translation;

        for (int i = 0; i < 1000; ++i) {
            parameter::BodyTransformParametersRelative zero_params({0, 0, 0}, {0, 0, 0});
            transformer->apply(zero_params, 0u);
        }

        auto& params = rigidbody.conformation->absolute_parameters.parameters[0];
        REQUIRE_THAT(params.rotation.x(), Catch::Matchers::WithinAbs(original_rotation.x(), 1e-10));
        REQUIRE_THAT(params.rotation.y(), Catch::Matchers::WithinAbs(original_rotation.y(), 1e-10));
        REQUIRE_THAT(params.rotation.z(), Catch::Matchers::WithinAbs(original_rotation.z(), 1e-10));
        REQUIRE_THAT(params.translation.x(), Catch::Matchers::WithinAbs(original_translation.x(), 1e-10));
        REQUIRE_THAT(params.translation.y(), Catch::Matchers::WithinAbs(original_translation.y(), 1e-10));
        REQUIRE_THAT(params.translation.z(), Catch::Matchers::WithinAbs(original_translation.z(), 1e-10));
    }

    SECTION("zero transforms after non-zero transform preserve state") {
        // First apply a non-trivial transform
        parameter::BodyTransformParametersRelative real_params({1.5, 2.3, -0.7}, {0.3, -0.1, 0.5});
        transformer->apply(real_params, 0u);

        auto rotation_after = rigidbody.conformation->absolute_parameters.parameters[0].rotation;
        auto translation_after = rigidbody.conformation->absolute_parameters.parameters[0].translation;

        // Then apply many zero transforms
        for (int i = 0; i < 1000; ++i) {
            parameter::BodyTransformParametersRelative zero_params({0, 0, 0}, {0, 0, 0});
            transformer->apply(zero_params, 0u);
        }

        auto& params = rigidbody.conformation->absolute_parameters.parameters[0];
        REQUIRE_THAT(params.rotation.x(), Catch::Matchers::WithinAbs(rotation_after.x(), 1e-10));
        REQUIRE_THAT(params.rotation.y(), Catch::Matchers::WithinAbs(rotation_after.y(), 1e-10));
        REQUIRE_THAT(params.rotation.z(), Catch::Matchers::WithinAbs(rotation_after.z(), 1e-10));
        REQUIRE_THAT(params.translation.x(), Catch::Matchers::WithinAbs(translation_after.x(), 1e-10));
        REQUIRE_THAT(params.translation.y(), Catch::Matchers::WithinAbs(translation_after.y(), 1e-10));
        REQUIRE_THAT(params.translation.z(), Catch::Matchers::WithinAbs(translation_after.z(), 1e-10));
    }
}

struct fixture {
    fixture() {
        settings::molecule::implicit_hydrogens = false;

        a1 = AtomFF({-1, -1, -1}, form_factor::form_factor_t::C);
        a2 = AtomFF({-1,  1, -1}, form_factor::form_factor_t::C);
        b1 = Body(std::vector{a1, a2});

        a3 = AtomFF({1, -1, -1}, form_factor::form_factor_t::C);
        a4 = AtomFF({1,  1, -1}, form_factor::form_factor_t::C);
        b2 = Body(std::vector{a3, a4});

        a5 =  AtomFF({-1, -1,  1}, form_factor::form_factor_t::C);
        a6 =  AtomFF({-1,  1,  1}, form_factor::form_factor_t::C);
        b3 = Body(std::vector{a5, a6});

        a7 =  AtomFF({ 1, -1,  1}, form_factor::form_factor_t::C);
        a8 =  AtomFF({ 1,  1,  1}, form_factor::form_factor_t::C);
        b4 = Body(std::vector{a7, a8});

        a9 =  AtomFF({ 0,  0,  0}, form_factor::form_factor_t::C);
        a10 = AtomFF({ 0,  0,  2}, form_factor::form_factor_t::C);
        b5 = Body(std::vector{a9, a10});

        bodies = {b1, b2, b3, b4, b5};
    }

    AtomFF a1, a2, a3, a4, a5, a6, a7, a8, a9, a10;
    Body b1, b2, b3, b4, b5;
    std::vector<Body> bodies;
};

TEST_CASE_METHOD(fixture, "TransformStrategy::apply") {
    settings::general::verbose = false;
    settings::grid::scaling = 2;

    SECTION("SingleTransform::apply") {
        Rigidbody rigidbody(Molecule{bodies}); test::mark_backbone_carbons(rigidbody.molecule);
        auto& manager = rigidbody.constraints;

        manager->add_constraint(std::make_unique<rigidbody::constraints::DistanceConstraintBond>(&rigidbody.molecule, 0, 1)); // 0 <-- this one
        manager->add_constraint(std::make_unique<rigidbody::constraints::DistanceConstraintBond>(&rigidbody.molecule, 1, 2)); // 1
        manager->add_constraint(std::make_unique<rigidbody::constraints::DistanceConstraintBond>(&rigidbody.molecule, 2, 3)); // 2
        manager->add_constraint(std::make_unique<rigidbody::constraints::DistanceConstraintBond>(&rigidbody.molecule, 3, 4)); // 3

        // translate
        rigidbody::transform::SingleTransform transform(&rigidbody);
        transform.apply({{1, 0, 0}, {0, 0, 0}}, manager->discoverable_constraints[0].get(), manager->discoverable_constraints[0]->ibody1);
        CHECK(rigidbody.molecule.get_body(0).get_atom(0).coordinates() == Vector3<double>(0, -1, -1));
        CHECK(rigidbody.molecule.get_body(0).get_atom(1).coordinates() == Vector3<double>(0,  1, -1));
        transform.undo();
        CHECK(rigidbody.molecule.get_body(0).get_atom(0).coordinates() == Vector3<double>(-1, -1, -1));
        CHECK(rigidbody.molecule.get_body(0).get_atom(1).coordinates() == Vector3<double>(-1,  1, -1));

        // rotate
        transform.apply({{0, 0, 0}, {0, 0, std::numbers::pi/2}}, manager->discoverable_constraints[0].get(), manager->discoverable_constraints[0]->ibody1);
        CHECK(rigidbody.molecule.get_body(0).get_atom(0).coordinates() == Vector3<double>( 1, -3, -1));
        CHECK(rigidbody.molecule.get_body(0).get_atom(1).coordinates() == Vector3<double>(-1, -3, -1));
        transform.undo();
        CHECK(rigidbody.molecule.get_body(0).get_atom(0).coordinates() == Vector3<double>(-1, -1, -1));
        CHECK(rigidbody.molecule.get_body(0).get_atom(1).coordinates() == Vector3<double>(-1,  1, -1));
    }

    SECTION("RigidTransform::apply") {
        // generate an easier constraint to test the apply function on
        AtomFF a11({3,  1, -1}, form_factor::form_factor_t::C);
        AtomFF a12({3, -1,  1}, form_factor::form_factor_t::C);
        Body b6(std::vector<AtomFF>{a11, a12});
        Rigidbody rigidbody(Molecule{std::vector<Body>{b1, b2, b3, b4, b5, b6}}); test::mark_backbone_carbons(rigidbody.molecule);
        auto& manager = rigidbody.constraints;

        manager->add_constraint(std::make_unique<rigidbody::constraints::DistanceConstraintBond>(&rigidbody.molecule, 0, 1)); // 0 <-- first this one
        manager->add_constraint(std::make_unique<rigidbody::constraints::DistanceConstraintBond>(&rigidbody.molecule, 1, 5)); // 1 <-- then this one
        manager->add_constraint(std::make_unique<rigidbody::constraints::DistanceConstraintBond>(&rigidbody.molecule, 5, 2)); // 2
        manager->add_constraint(std::make_unique<rigidbody::constraints::DistanceConstraintBond>(&rigidbody.molecule, 2, 3)); // 3
        manager->add_constraint(std::make_unique<rigidbody::constraints::DistanceConstraintBond>(&rigidbody.molecule, 3, 4)); // 4

        // single-body translate
        rigidbody::transform::RigidTransform transform(&rigidbody);
        transform.apply({{1, 0, 0}, {0, 0, 0}}, manager->discoverable_constraints[0].get(), manager->discoverable_constraints[0]->ibody1);
        CHECK(rigidbody.molecule.get_body(0).get_atom(0).coordinates() == Vector3<double>(0, -1, -1));
        CHECK(rigidbody.molecule.get_body(0).get_atom(1).coordinates() == Vector3<double>(0,  1, -1));
        transform.undo();
        CHECK(rigidbody.molecule.get_body(0).get_atom(0).coordinates() == Vector3<double>(-1, -1, -1));
        CHECK(rigidbody.molecule.get_body(0).get_atom(1).coordinates() == Vector3<double>(-1,  1, -1));

        // single-body rotate
        transform.apply({{0, 0, 0}, {0, 0, std::numbers::pi/2}}, manager->discoverable_constraints[0].get(), manager->discoverable_constraints[0]->ibody1);
        CHECK(rigidbody.molecule.get_body(0).get_atom(0).coordinates() == Vector3<double>( 1, -3, -1));
        CHECK(rigidbody.molecule.get_body(0).get_atom(1).coordinates() == Vector3<double>(-1, -3, -1));
        transform.undo();
        CHECK(rigidbody.molecule.get_body(0).get_atom(0).coordinates() == Vector3<double>(-1, -1, -1));
        CHECK(rigidbody.molecule.get_body(0).get_atom(1).coordinates() == Vector3<double>(-1,  1, -1));

        // multi-body translate
        transform.apply({{1, 0, 0}, {0, 0, 0}}, manager->discoverable_constraints[1].get(), manager->discoverable_constraints[1]->ibody1);
        CHECK(rigidbody.molecule.get_body(0).get_atom(0).coordinates() == Vector3<double>(0, -1, -1));
        CHECK(rigidbody.molecule.get_body(0).get_atom(1).coordinates() == Vector3<double>(0,  1, -1));
        CHECK(rigidbody.molecule.get_body(1).get_atom(0).coordinates() == Vector3<double>(2, -1, -1));
        CHECK(rigidbody.molecule.get_body(1).get_atom(1).coordinates() == Vector3<double>(2,  1, -1));
        transform.undo();
        CHECK(rigidbody.molecule.get_body(0).get_atom(0).coordinates() == Vector3<double>(-1, -1, -1));
        CHECK(rigidbody.molecule.get_body(0).get_atom(1).coordinates() == Vector3<double>(-1,  1, -1));
        CHECK(rigidbody.molecule.get_body(1).get_atom(0).coordinates() == Vector3<double>( 1, -1, -1));
        CHECK(rigidbody.molecule.get_body(1).get_atom(1).coordinates() == Vector3<double>( 1,  1, -1));

        // multi-body rotate
        transform.apply({{0, 0, 0}, {0, 0, std::numbers::pi/2}}, manager->discoverable_constraints[1].get(), manager->discoverable_constraints[1]->ibody1);
        CHECK(rigidbody.molecule.get_body(0).get_atom(0).coordinates() == Vector3<double>( 5, -3, -1));
        CHECK(rigidbody.molecule.get_body(0).get_atom(1).coordinates() == Vector3<double>( 3, -3, -1));
        CHECK(rigidbody.molecule.get_body(1).get_atom(0).coordinates() == Vector3<double>( 5, -1, -1));
        CHECK(rigidbody.molecule.get_body(1).get_atom(1).coordinates() == Vector3<double>( 3, -1, -1));
        transform.undo();
        CHECK(rigidbody.molecule.get_body(0).get_atom(0).coordinates() == Vector3<double>(-1, -1, -1));
        CHECK(rigidbody.molecule.get_body(0).get_atom(1).coordinates() == Vector3<double>(-1,  1, -1));
        CHECK(rigidbody.molecule.get_body(1).get_atom(0).coordinates() == Vector3<double>( 1, -1, -1));
        CHECK(rigidbody.molecule.get_body(1).get_atom(1).coordinates() == Vector3<double>( 1,  1, -1));
    }

    SECTION("Parameters can reconstruct body after transformation") {
        Rigidbody rigidbody(Molecule{bodies}); test::mark_backbone_carbons(rigidbody.molecule);
        auto& manager = rigidbody.constraints;

        manager->add_constraint(std::make_unique<rigidbody::constraints::DistanceConstraintBond>(&rigidbody.molecule, 0, 1));
        manager->add_constraint(std::make_unique<rigidbody::constraints::DistanceConstraintBond>(&rigidbody.molecule, 1, 2));
        manager->add_constraint(std::make_unique<rigidbody::constraints::DistanceConstraintBond>(&rigidbody.molecule, 2, 3));
        manager->add_constraint(std::make_unique<rigidbody::constraints::DistanceConstraintBond>(&rigidbody.molecule, 3, 4));

        rigidbody::transform::RigidTransform transform(&rigidbody);

        // Apply a rotation then check that parameters can reconstruct
        transform.apply({{0, 0, 0}, {0, 0, std::numbers::pi/2}}, manager->discoverable_constraints[0].get(), manager->discoverable_constraints[0]->ibody1);

        for (int ibody = 0; ibody < rigidbody.molecule.size_body(); ++ibody) {
            auto& current_body = rigidbody.molecule.get_body(ibody);
            auto& params = rigidbody.conformation->absolute_parameters.parameters[ibody];
            auto& original = rigidbody.conformation->initial_conformation[ibody];

            // Reconstruct from original + parameters
            Body reconstructed = original;
            reconstructed.rotate(matrix::rotation_matrix(params.rotation));
            reconstructed.translate(params.translation);

            auto current_cm = current_body.get_cm();
            auto reconstructed_cm = reconstructed.get_cm();

            INFO("Body " << ibody << " should be reconstructible from parameters");
            CHECK(std::abs(reconstructed_cm.x() - current_cm.x()) < 1e-6);
            CHECK(std::abs(reconstructed_cm.y() - current_cm.y()) < 1e-6);
            CHECK(std::abs(reconstructed_cm.z() - current_cm.z()) < 1e-6);
        }
    }
}

static auto vector_contains = [] (std::vector<int> vec, const std::vector<int>& vals) {
    std::unordered_set<int> set(vec.begin(), vec.end());
    return std::ranges::all_of(vals, [&set] (int val) {return set.contains(val);});
};

TEST_CASE_METHOD(fixture, "RigidTransform::get_connected") {
    settings::general::verbose = false;

    SECTION("get_connected") {
        struct TestRigidTransform : public rigidbody::transform::RigidTransform {
            using RigidTransform::RigidTransform;
            using RigidTransform::get_connected;
        };

        SECTION("simple") {
            Rigidbody rigidbody(Molecule{bodies}); test::mark_backbone_carbons(rigidbody.molecule);
            auto& manager = rigidbody.constraints;

            manager->add_constraint(std::make_unique<rigidbody::constraints::DistanceConstraintBond>(&rigidbody.molecule, 0, 1)); // 0
            manager->add_constraint(std::make_unique<rigidbody::constraints::DistanceConstraintBond>(&rigidbody.molecule, 1, 2)); // 1
            manager->add_constraint(std::make_unique<rigidbody::constraints::DistanceConstraintBond>(&rigidbody.molecule, 2, 3)); // 2
            manager->add_constraint(std::make_unique<rigidbody::constraints::DistanceConstraintBond>(&rigidbody.molecule, 3, 4)); // 3

            TestRigidTransform transform(&rigidbody);

            // 0 - 1 - 2 - 3 - 4            //
            auto group1 = transform.get_connected(rigidbody.constraints->discoverable_constraints[0].get());
            REQUIRE(group1.indices.size() == 1);
            CHECK(group1.indices[0] == 0);

            auto group2 = transform.get_connected(rigidbody.constraints->discoverable_constraints[1].get());
            REQUIRE(group2.indices.size() == 2);
            CHECK(vector_contains(group2.indices, {0, 1}));

            auto group3 = transform.get_connected(rigidbody.constraints->discoverable_constraints[2].get());
            REQUIRE(group3.indices.size() == 2);
            CHECK(vector_contains(group3.indices, {3, 4}));

            auto group4 = transform.get_connected(rigidbody.constraints->discoverable_constraints[3].get());
            REQUIRE(group4.indices.size() == 1);
            CHECK(group4.indices[0] == 4);
        }

        SECTION("complex") {
            Body b6(std::vector<AtomFF>{AtomFF({0, 1, 0}, form_factor::form_factor_t::C), AtomFF({0, 2, 0}, form_factor::form_factor_t::C)});
            Body b7(std::vector<AtomFF>{AtomFF({0, 3, 0}, form_factor::form_factor_t::C), AtomFF({0, 4, 0}, form_factor::form_factor_t::C)});
            Body b8(std::vector<AtomFF>{AtomFF({1, 0, 0}, form_factor::form_factor_t::C), AtomFF({2, 0, 0}, form_factor::form_factor_t::C)});
            bodies = {b1, b2, b3, b4, b5, b6, b7, b8};
            Rigidbody rigidbody(Molecule{bodies}); test::mark_backbone_carbons(rigidbody.molecule);
            auto& manager = rigidbody.constraints;

            manager->add_constraint(std::make_unique<rigidbody::constraints::DistanceConstraintBond>(&rigidbody.molecule, 0, 1)); // 0
            manager->add_constraint(std::make_unique<rigidbody::constraints::DistanceConstraintBond>(&rigidbody.molecule, 1, 2)); // 1
            manager->add_constraint(std::make_unique<rigidbody::constraints::DistanceConstraintBond>(&rigidbody.molecule, 2, 3)); // 2
            manager->add_constraint(std::make_unique<rigidbody::constraints::DistanceConstraintBond>(&rigidbody.molecule, 3, 4)); // 3
            manager->add_constraint(std::make_unique<rigidbody::constraints::DistanceConstraintBond>(&rigidbody.molecule, 3, 5)); // 4
            manager->add_constraint(std::make_unique<rigidbody::constraints::DistanceConstraintBond>(&rigidbody.molecule, 5, 6)); // 5
            manager->add_constraint(std::make_unique<rigidbody::constraints::DistanceConstraintBond>(&rigidbody.molecule, 3, 7)); // 6

            //             5 - 6            //
            //             |                //
            // 0 - 1 - 2 - 3 - 4            //
            //             |                //
            //             7                //

            TestRigidTransform transform(&rigidbody);
            auto group1 = transform.get_connected(rigidbody.constraints->discoverable_constraints[0].get());
            REQUIRE(group1.indices.size() == 1);
            CHECK(group1.indices[0] == 0);

            auto group2 = transform.get_connected(rigidbody.constraints->discoverable_constraints[1].get());
            REQUIRE(group2.indices.size() == 2);
            CHECK(vector_contains(group2.indices, {0, 1}));

            auto group3 = transform.get_connected(rigidbody.constraints->discoverable_constraints[2].get());
            REQUIRE(group3.indices.size() == 3);
            CHECK(vector_contains(group3.indices, {0, 1, 2}));

            auto group4 = transform.get_connected(rigidbody.constraints->discoverable_constraints[3].get());
            REQUIRE(group4.indices.size() == 1);
            CHECK(group4.indices[0] == 4);

            auto group5 = transform.get_connected(rigidbody.constraints->discoverable_constraints[4].get());
            REQUIRE(group5.indices.size() == 2);
            CHECK(vector_contains(group5.indices, {5, 6}));

            auto group6 = transform.get_connected(rigidbody.constraints->discoverable_constraints[5].get());
            REQUIRE(group6.indices.size() == 1);
            CHECK(group6.indices[0] == 6);

            auto group7 = transform.get_connected(rigidbody.constraints->discoverable_constraints[6].get());
            REQUIRE(group7.indices.size() == 1);
            CHECK(group7.indices[0] == 7);
        }
    }
}

TEST_CASE_METHOD(fixture, "SingleTransform::apply: moves the selected body, whichever end of the constraint it is") {
    settings::general::verbose = false;
    settings::grid::scaling = 2;

    // declared "b2-b1": body 1 is the constraint's ibody1, so HEAD moved body 1 whichever body was selected
    Rigidbody rigidbody(Molecule{bodies}); test::mark_backbone_carbons(rigidbody.molecule);
    auto& manager = rigidbody.constraints;
    manager->add_constraint(std::make_unique<rigidbody::constraints::DistanceConstraintBond>(&rigidbody.molecule, 1, 0));
    auto* constraint = manager->discoverable_constraints[0].get();
    REQUIRE(constraint->ibody1 == 1);

    auto coords = [&] (int ibody) {
        std::vector<Vector3<double>> out;
        for (const auto& a : rigidbody.molecule.get_body(ibody).get_atoms()) {out.push_back(a.coordinates());}
        return out;
    };
    rigidbody::transform::SingleTransform transform(&rigidbody);

    for (int selected : {0, 1}) {
        int other = 1 - selected;
        auto before_sel = coords(selected), before_other = coords(other);
        Vector3<double> pivot = (constraint->ibody1 == selected ? constraint->get_atom2() : constraint->get_atom1()).coordinates();

        // translate: the selected body moves, the other one does not
        transform.apply({{1, 0, 0}, {0, 0, 0}}, constraint, selected);
        auto after = coords(selected);
        for (std::size_t i = 0; i < after.size(); ++i) {CHECK(after[i] == before_sel[i] + Vector3<double>(1, 0, 0));}
        CHECK(coords(other) == before_other);
        transform.undo();
        CHECK(coords(selected) == before_sel);

        // rotate: about the constraint's atom on the other body, so every distance to it is kept
        transform.apply({{0, 0, 0}, {0, 0, std::numbers::pi/2}}, constraint, selected);
        after = coords(selected);
        CHECK(after != before_sel);
        for (std::size_t i = 0; i < after.size(); ++i) {
            CHECK_THAT((after[i] - pivot).magnitude(), Catch::Matchers::WithinAbs((before_sel[i] - pivot).magnitude(), 1e-9));
        }
        CHECK(coords(other) == before_other);
        transform.undo();
        CHECK(coords(selected) == before_sel);
    }
}
