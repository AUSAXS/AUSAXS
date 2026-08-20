#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rigidbody/Rigidbody.h>
#include <rigidbody/transform/TransformStrategy.h>
#include <rigidbody/transform/TransformGroup.h>
#include <rigidbody/constraints/DistanceConstraintBond.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <data/Body.h>
#include <data/Molecule.h>
#include <data/detail/SimpleBody.h>
#include <rigidbody/parameters/BodyTransformParametersRelative.h>
#include <data/symmetry/PointSymmetry.h>
#include <data/symmetry/CyclicSymmetry.h>
#include <math/MatrixUtils.h>
#include <settings/All.h>

#include <numbers>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;

TEST_CASE("TransformStrategy::apply unconstrained body") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
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
    settings::molecule::center = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
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
    settings::molecule::center = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
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
    settings::molecule::center = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
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
        for (size_t i = 0; i < current_body.size_atom(); ++i) {
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
    settings::molecule::center = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
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
    settings::molecule::center = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
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
        std::sort(distances.begin(), distances.end());
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
