// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/elements/setup/ConvertToSymmetryElement.h>
#include <rigidbody/Rigidbody.h>
#include <data/symmetry/CyclicSymmetry.h>
#include <data/symmetry/CompositeSymmetry.h>
#include <data/symmetry/PointSymmetry.h>
#include <data/symmetry/PredefinedSymmetries.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <math/MatrixUtils.h>
#include <settings/All.h>
#include <io/ExistingFile.h>

#include <cstdio>
#include <fstream>
#include <numbers>

using namespace ausaxs;
using namespace ausaxs::symmetry;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::sequencer;

namespace {
    struct Fixture {
        Fixture() {
            settings::general::verbose = false;
            settings::molecule::implicit_hydrogens = false;
            settings::grid::min_bins = 20;
        }
    };

    // reference atom cloud (non-degenerate)
    std::vector<Vector3<double>> reference_atoms() {
        return {{2, 0, 0}, {3, 1, 0.5}, {2.5, -1, 1}, {4, 0.5, -0.5}, {2, 2, 1}, {5, -1, 0}};
    }

    void write_pdb(const std::string& path, const std::vector<Vector3<double>>& atoms) {
        std::ofstream out(path);
        int serial = 1;
        for (const auto& a : atoms) {
            char line[128];
            std::snprintf(line, sizeof line,
                "ATOM  %5d  C   ALA A%4d    %8.3f%8.3f%8.3f  1.00  0.00           C\n",
                serial, serial, a.x(), a.y(), a.z());
            out << line;
            ++serial;
        }
        out << "END\n";
    }
}

TEST_CASE_METHOD(Fixture, "ConvertToSymmetryElement collapses a cyclic assembly") {
    // three bodies related by a C3 rotation about an arbitrary centre/axis
    auto ref = reference_atoms();
    CyclicSymmetry source({{1, 2, 0}}, {{0, 0, 0}}, {0.2, 1, 0.4}, 2*std::numbers::pi/3, 2);
    Vector3<double> cm{0, 0, 0};
    for (const auto& a : ref) {cm += a;}
    cm /= static_cast<double>(ref.size());

    std::vector<std::vector<Vector3<double>>> chains{ref};
    for (int k = 1; k <= 2; ++k) {
        auto t = source.get_transform(cm, k);
        std::vector<Vector3<double>> chain;
        for (const auto& a : ref) {chain.push_back(t(a));}
        chains.push_back(std::move(chain));
    }

    std::vector<std::string> files = {
        "/tmp/ausaxs_convsym_0.pdb", "/tmp/ausaxs_convsym_1.pdb", "/tmp/ausaxs_convsym_2.pdb"
    };
    for (std::size_t i = 0; i < files.size(); ++i) {write_pdb(files[i], chains[i]);}

    Sequencer seq(io::ExistingFile("tests/files/SASDJG5.dat"));
    seq.setup().load(files);
    auto* molecule = seq._get_molecule();
    REQUIRE(molecule->size_body() == 3);

    // collapse the three bodies into body 0 + a fitted C3
    ConvertToSymmetryElement convert(&seq, {0, 1, 2}, "c3");

    // only the reference body remains, now carrying one symmetry generating two copies
    REQUIRE(molecule->size_body() == 1);
    REQUIRE(molecule->get_body(0).size_symmetry() == 1);
    CHECK(molecule->get_body(0).symmetry().get(0)->repetitions() == 2);

    // expanding the collapsed body reproduces the original three-chain assembly
    auto expanded = molecule->get_body(0).symmetry().explicit_structure();
    REQUIRE(expanded.atoms.size() == 3*ref.size());
    for (std::size_t copy = 0; copy < 3; ++copy) {
        for (std::size_t i = 0; i < ref.size(); ++i) {
            Vector3<double> got = expanded.atoms[copy*ref.size() + i].coordinates();
            // tolerance reflects the PDB writer's 3-decimal coordinate rounding, not fit accuracy
            CHECK_THAT((got - chains[copy][i]).magnitude(), Catch::Matchers::WithinAbs(0, 5e-3));
        }
    }
}

TEST_CASE_METHOD(Fixture, "ConvertToSymmetryElement collapses a composite p2-p2 assembly") {
    // four bodies from an inner p2 dimer replicated by an outer p2 (the A2M-native use case)
    auto ref = reference_atoms();
    auto source = symmetry::create("p2-p2");
    auto* comp = dynamic_cast<symmetry::CompositeSymmetry*>(source.get());
    REQUIRE(comp != nullptr);
    dynamic_cast<symmetry::PointSymmetry*>(comp->inner.get())->translation = {5, 1, -2};
    dynamic_cast<symmetry::PointSymmetry*>(comp->inner.get())->rotation = {0.3, 0.7, -0.4};
    dynamic_cast<symmetry::PointSymmetry*>(comp->outer.get())->translation = {-3, 6, 1};
    dynamic_cast<symmetry::PointSymmetry*>(comp->outer.get())->rotation = {1.1, -0.2, 0.5};

    Vector3<double> cm{0, 0, 0};
    for (const auto& a : ref) {cm += a;}
    cm /= static_cast<double>(ref.size());

    std::vector<std::vector<Vector3<double>>> chains{ref};
    for (int k = 1; k <= 3; ++k) {
        auto t = source->get_transform(cm, k);
        std::vector<Vector3<double>> chain;
        for (const auto& a : ref) {chain.push_back(t(a));}
        chains.push_back(std::move(chain));
    }

    std::vector<std::string> files;
    for (int i = 0; i < 4; ++i) {
        files.push_back("/tmp/ausaxs_convp2_" + std::to_string(i) + ".pdb");
        write_pdb(files.back(), chains[i]);
    }

    Sequencer seq(io::ExistingFile("tests/files/SASDJG5.dat"));
    seq.setup().load(files);
    auto* molecule = seq._get_molecule();
    REQUIRE(molecule->size_body() == 4);

    ConvertToSymmetryElement convert(&seq, {0, 1, 2, 3}, "p2-p2");

    REQUIRE(molecule->size_body() == 1);
    REQUIRE(molecule->get_body(0).size_symmetry() == 1);
    CHECK(molecule->get_body(0).symmetry().get(0)->repetitions() == 3);

    auto expanded = molecule->get_body(0).symmetry().explicit_structure();
    REQUIRE(expanded.atoms.size() == 4*ref.size());
    for (std::size_t copy = 0; copy < 4; ++copy) {
        for (std::size_t i = 0; i < ref.size(); ++i) {
            Vector3<double> got = expanded.atoms[copy*ref.size() + i].coordinates();
            CHECK_THAT((got - chains[copy][i]).magnitude(), Catch::Matchers::WithinAbs(0, 5e-3));
        }
    }
}

TEST_CASE_METHOD(Fixture, "ConvertToSymmetryElement rejects an assembly that is not symmetric") {
    // three bodies that are pure translations of each other do not form a c3, so the fit must fail
    auto ref = reference_atoms();
    std::vector<std::vector<Vector3<double>>> chains{ref};
    chains.push_back([&]{ auto c = ref; for (auto& p : c) {p += Vector3<double>{15, 0, 0};} return c; }());
    chains.push_back([&]{ auto c = ref; for (auto& p : c) {p += Vector3<double>{0, 15, 0};} return c; }());

    std::vector<std::string> files;
    for (int i = 0; i < 3; ++i) {
        files.push_back("/tmp/ausaxs_convbad_" + std::to_string(i) + ".pdb");
        write_pdb(files.back(), chains[i]);
    }

    Sequencer seq(io::ExistingFile("tests/files/SASDJG5.dat"));
    seq.setup().load(files);
    CHECK_THROWS([&]{ ConvertToSymmetryElement convert(&seq, {0, 1, 2}, "c3"); }());
}

TEST_CASE_METHOD(Fixture, "ConvertToSymmetryElement is reachable through the script parser") {
    // the inline "convert_to_symmetry c3" form, driven end-to-end through the parser
    auto ref = reference_atoms();
    CyclicSymmetry source({{1, 2, 0}}, {{0, 0, 0}}, {0.2, 1, 0.4}, 2*std::numbers::pi/3, 2);
    Vector3<double> cm{0, 0, 0};
    for (const auto& a : ref) {cm += a;}
    cm /= static_cast<double>(ref.size());

    std::vector<std::string> files;
    for (int k = 0; k < 3; ++k) {
        std::vector<Vector3<double>> chain;
        if (k == 0) {chain = ref;}
        else {auto t = source.get_transform(cm, k); for (const auto& a : ref) {chain.push_back(t(a));}}
        files.push_back("/tmp/ausaxs_convparse_" + std::to_string(k) + ".pdb");
        write_pdb(files.back(), chain);
    }

    std::string script =
        "load {\n"
        "pdb " + files[0] + " " + files[1] + " " + files[2] + "\n"
        "saxs tests/files/SASDJG5.dat\n"
        "}\n"
        "convert_to_symmetry c3\n";

    auto seq = SequenceParser().parse_text(script);
    auto* molecule = seq->_get_molecule();
    REQUIRE(molecule->size_body() == 1);
    REQUIRE(molecule->get_body(0).size_symmetry() == 1);
    CHECK(molecule->get_body(0).symmetry().get(0)->repetitions() == 2);
}

TEST_CASE_METHOD(Fixture, "ConvertToSymmetryElement on a real p2 dimer (SASDJG5)") {
    // SASDJG5 is a two-chain p2 homodimer
    Sequencer seq(io::ExistingFile("tests/files/SASDJG5.dat"));
    seq.setup().load("tests/files/SASDJG5.pdb");
    auto* molecule = seq._get_molecule();
    REQUIRE(molecule->size_body() == 2);

    REQUIRE_NOTHROW([&]{ ConvertToSymmetryElement convert(&seq, {0, 1}, "p2"); }());
    REQUIRE(molecule->size_body() == 1);
    REQUIRE(molecule->get_body(0).size_symmetry() == 1);
    CHECK(molecule->get_body(0).symmetry().get(0)->repetitions() == 1);
}

TEST_CASE_METHOD(Fixture, "ConvertToSymmetryElement on a real p2-p2 assembly (A2M_native)", "[slow]") {
    // A2M-native is a four-chain assembly with p2-p2 symmetry; chains may be in any order, so this
    // also exercises the ordering search on real data
    Sequencer seq(io::ExistingFile("tests/files/A2M_native.dat"));
    seq.setup().load("tests/files/A2M_native.pdb");
    auto* molecule = seq._get_molecule();
    REQUIRE(molecule->size_body() == 4);

    REQUIRE_NOTHROW([&]{ ConvertToSymmetryElement convert(&seq, {0, 1, 2, 3}, "p2-p2"); }());
    REQUIRE(molecule->size_body() == 1);
    REQUIRE(molecule->get_body(0).size_symmetry() == 1);
    CHECK(molecule->get_body(0).symmetry().get(0)->repetitions() == 3);
}

TEST_CASE_METHOD(Fixture, "ConvertToSymmetryElement rejects a wrong body count") {
    auto ref = reference_atoms();
    std::vector<std::string> files = {"/tmp/ausaxs_convsym_a.pdb", "/tmp/ausaxs_convsym_b.pdb"};
    write_pdb(files[0], ref);
    write_pdb(files[1], ref);

    Sequencer seq(io::ExistingFile("tests/files/SASDJG5.dat"));
    seq.setup().load(files);

    // c3 needs three bodies, only two are available
    CHECK_THROWS([&]{ ConvertToSymmetryElement convert(&seq, {0, 1}, "c3"); }());
}
