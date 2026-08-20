// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/elements/setup/ConvertToSymmetryElement.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/detail/SystemSpecification.h>
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
#include <set>

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

    // The reference cloud plus its two images under a C3 about an arbitrary centre/axis, reference first.
    std::vector<std::vector<Vector3<double>>> c3_assembly() {
        auto ref = reference_atoms();
        CyclicSymmetry source({{1, 2, 0}}, {{0, 0, 0}}, {0.2, 1, 0.4}, 2*std::numbers::pi/3, 2);
        Vector3<double> cm{0, 0, 0};
        for (const auto& a : ref) {cm += a;}
        cm /= static_cast<double>(ref.size());

        std::vector<std::vector<Vector3<double>>> chains{ref};
        for (int k = 1; k <= 2; ++k) {
            auto t = source._get_transform(cm, k);
            std::vector<Vector3<double>> chain;
            for (const auto& a : ref) {chain.push_back(t(a));}
            chains.push_back(std::move(chain));
        }
        return chains;
    }

    // Write an atom cloud as a single chain. Each atom is its own residue unless `residues` assigns them explicitly, and the atom indices in `omit` are left
    // out entirely - which, since the residue ids stay tied to the atom's index in the reference cloud, is how a chain modelled to a lesser extent than its
    // copies appears in a deposited structure: the retained atoms keep the ids they would have had.
    void write_pdb(
        const io::File& path, const std::vector<Vector3<double>>& atoms,
        const std::vector<int>& residues = {}, const std::set<std::size_t>& omit = {}
    ) {
        path.create();
        std::ofstream out(path);
        int serial = 1;
        for (std::size_t i = 0; i < atoms.size(); ++i) {
            if (omit.contains(i)) {continue;}
            char line[128];
            std::snprintf(line, sizeof line,
                "ATOM  %5d  C   ALA A%4d    %8.3f%8.3f%8.3f  1.00  0.00           C\n",
                serial, residues.empty() ? static_cast<int>(i) + 1 : residues[i], atoms[i].x(), atoms[i].y(), atoms[i].z());
            out << line;
            ++serial;
        }
        out << "END\n";
    }
}

TEST_CASE_METHOD(Fixture, "ConvertToSymmetryElement collapses a cyclic assembly") {
    // three bodies related by a C3 rotation about an arbitrary centre/axis
    auto ref = reference_atoms();
    auto chains = c3_assembly();

    std::vector<std::string> files = {
        "temp/rigidbody/ausaxs_convsym_0.pdb", "temp/rigidbody/ausaxs_convsym_1.pdb", "temp/rigidbody/ausaxs_convsym_2.pdb"
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
        auto t = source->_get_transform(cm, k);
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

TEST_CASE_METHOD(Fixture, "ConvertToSymmetryElement matches up copies modelled to differing extents") {
    // The copies are missing residue 3, in the middle of the chain: the correspondence must be recovered from the residue ids, since pairing the atoms off by
    // their position in the atom vector would misalign everything past the gap. A trailing gap - the more common deposited case - is the easier half of this.
    auto chains = c3_assembly();
    const auto& ref = chains[0];

    std::vector<std::string> files;
    for (int i = 0; i < 3; ++i) {
        files.push_back("temp/rigidbody/ausaxs_convgap_" + std::to_string(i) + ".pdb");
        write_pdb(files.back(), chains[i], {}, i == 0 ? std::set<std::size_t>{} : std::set<std::size_t>{2});
    }

    Sequencer seq(io::ExistingFile("tests/files/SASDJG5.dat"));
    seq.setup().load(files);
    auto* molecule = seq._get_molecule();
    REQUIRE(molecule->size_body() == 3);
    REQUIRE(molecule->get_body(0).size_atom() == ref.size());
    REQUIRE(molecule->get_body(1).size_atom() == ref.size() - 1);

    ConvertToSymmetryElement convert(&seq, {0, 1, 2}, "c3");

    // the fit only saw the five shared residues, but the primary body is kept whole and is what the symmetry replicates
    REQUIRE(molecule->size_body() == 1);
    REQUIRE(molecule->get_body(0).size_atom() == ref.size());
    REQUIRE(molecule->get_body(0).size_symmetry() == 1);

    auto expanded = molecule->get_body(0).symmetry().explicit_structure();
    REQUIRE(expanded.atoms.size() == 3*ref.size());
    for (std::size_t copy = 0; copy < 3; ++copy) {
        for (std::size_t i = 0; i < ref.size(); ++i) {
            // every atom is checked, including the one the copies never modelled: the symmetry regenerates it from the reference
            Vector3<double> got = expanded.atoms[copy*ref.size() + i].coordinates();
            CHECK_THAT((got - chains[copy][i]).magnitude(), Catch::Matchers::WithinAbs(0, 5e-3));
        }
    }
}

TEST_CASE_METHOD(Fixture, "ConvertToSymmetryElement drops residues modelled to differing extents") {
    // Three residues of two atoms each, with one copy missing the second half of residue 2. That residue offers no correspondence in either direction, so it
    // must be dropped from the fit entirely rather than matched atom-for-atom against the reference's two.
    auto chains = c3_assembly();
    const auto& ref = chains[0];
    std::vector<int> residues{1, 1, 2, 2, 3, 3};

    std::vector<std::string> files;
    for (int i = 0; i < 3; ++i) {
        files.push_back("temp/rigidbody/ausaxs_convpartial_" + std::to_string(i) + ".pdb");
        write_pdb(files.back(), chains[i], residues, i == 1 ? std::set<std::size_t>{3} : std::set<std::size_t>{});
    }

    Sequencer seq(io::ExistingFile("tests/files/SASDJG5.dat"));
    seq.setup().load(files);
    auto* molecule = seq._get_molecule();
    REQUIRE(molecule->size_body() == 3);

    ConvertToSymmetryElement convert(&seq, {0, 1, 2}, "c3");

    REQUIRE(molecule->size_body() == 1);
    REQUIRE(molecule->get_body(0).size_atom() == ref.size());
    auto expanded = molecule->get_body(0).symmetry().explicit_structure();
    REQUIRE(expanded.atoms.size() == 3*ref.size());
    for (std::size_t copy = 0; copy < 3; ++copy) {
        for (std::size_t i = 0; i < ref.size(); ++i) {
            Vector3<double> got = expanded.atoms[copy*ref.size() + i].coordinates();
            CHECK_THAT((got - chains[copy][i]).magnitude(), Catch::Matchers::WithinAbs(0, 5e-3));
        }
    }
}

TEST_CASE_METHOD(Fixture, "ConvertToSymmetryElement rejects copies that share no residues") {
    // a perfect c3 assembly, but with residue numbering that does not agree between the copies, leaving no correspondence to fit on
    auto chains = c3_assembly();

    std::vector<std::string> files;
    for (int i = 0; i < 3; ++i) {
        files.push_back("temp/rigidbody/ausaxs_convdisjoint_" + std::to_string(i) + ".pdb");
        std::vector<int> residues;
        for (std::size_t j = 0; j < chains[i].size(); ++j) {residues.push_back(100*i + static_cast<int>(j) + 1);}
        write_pdb(files.back(), chains[i], residues);
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
        else {auto t = source._get_transform(cm, k); for (const auto& a : ref) {chain.push_back(t(a));}}
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

TEST_CASE_METHOD(Fixture, "ConvertToSymmetryElement splits a single assembled body") {
    // the same C3 assembly as above, but written as a *single* file: the three copies share one body, so the element must decompose it itself
    auto ref = reference_atoms();
    CyclicSymmetry source({{1, 2, 0}}, {{0, 0, 0}}, {0.2, 1, 0.4}, 2*std::numbers::pi/3, 2);
    Vector3<double> cm{0, 0, 0};
    for (const auto& a : ref) {cm += a;}
    cm /= static_cast<double>(ref.size());

    std::vector<std::vector<Vector3<double>>> chains{ref};
    std::vector<Vector3<double>> assembly = ref;
    for (int k = 1; k <= 2; ++k) {
        auto t = source._get_transform(cm, k);
        std::vector<Vector3<double>> chain;
        for (const auto& a : ref) {chain.push_back(t(a));}
        assembly.insert(assembly.end(), chain.begin(), chain.end());
        chains.push_back(std::move(chain));
    }

    io::File file("temp/rigidbody/ausaxs_convsym_single.pdb");
    write_pdb(file, assembly);

    Sequencer seq(io::ExistingFile("tests/files/SASDJG5.dat"));
    seq.setup().load(std::vector<std::string>{file});
    auto* molecule = seq._get_molecule();
    REQUIRE(molecule->size_body() == 1);
    REQUIRE(molecule->get_body(0).size_atom() == 3*ref.size());

    ConvertToSymmetryElement convert(&seq, {0}, "c3");

    // the body has been reduced to the first copy, and carries the symmetry generating the other two
    REQUIRE(molecule->size_body() == 1);
    CHECK(molecule->get_body(0).size_atom() == ref.size());
    REQUIRE(molecule->get_body(0).size_symmetry() == 1);
    CHECK(molecule->get_body(0).symmetry().get(0)->repetitions() == 2);

    // expanding it reproduces the original assembly
    auto expanded = molecule->get_body(0).symmetry().explicit_structure();
    REQUIRE(expanded.atoms.size() == 3*ref.size());
    for (std::size_t copy = 0; copy < 3; ++copy) {
        for (std::size_t i = 0; i < ref.size(); ++i) {
            Vector3<double> got = expanded.atoms[copy*ref.size() + i].coordinates();
            CHECK_THAT((got - chains[copy][i]).magnitude(), Catch::Matchers::WithinAbs(0, 5e-3));
        }
    }

    // the stored initial conformation must stay parallel-indexed to the live body and origin-centred, with the translation restoring the live position
    const auto& initial = seq._get_rigidbody()->conformation->initial_conformation[0];
    REQUIRE(initial.size_atom() == ref.size());
    CHECK_THAT(initial.get_cm().magnitude(), Catch::Matchers::WithinAbs(0, 1e-9));
    auto translation = seq._get_rigidbody()->conformation->absolute_parameters.parameters[0].translation;
    CHECK_THAT((translation - molecule->get_body(0).get_cm()).magnitude(), Catch::Matchers::WithinAbs(0, 1e-9));
}

TEST_CASE_METHOD(Fixture, "ConvertToSymmetryElement rejects a single body that does not divide evenly") {
    // 7 atoms cannot be split into the 3 copies a c3 needs
    std::vector<Vector3<double>> atoms;
    for (int i = 0; i < 7; ++i) {atoms.push_back({static_cast<double>(i), 0.5*i, static_cast<double>(-i)});}

    io::File file("temp/rigidbody/ausaxs_convsym_indivisible.pdb");
    write_pdb(file, atoms);

    Sequencer seq(io::ExistingFile("tests/files/SASDJG5.dat"));
    seq.setup().load(std::vector<std::string>{file});
    REQUIRE(seq._get_molecule()->size_body() == 1);

    CHECK_THROWS([&]{ ConvertToSymmetryElement convert(&seq, {0}, "c3"); }());
}

TEST_CASE_METHOD(Fixture, "ConvertToSymmetryElement rejects a single body that is not symmetric") {
    // three copies related by pure translation do not form a c3, so the split must not rescue them
    auto ref = reference_atoms();
    std::vector<Vector3<double>> assembly = ref;
    for (const auto& shift : {Vector3<double>{15, 0, 0}, Vector3<double>{0, 15, 0}}) {
        for (const auto& a : ref) {assembly.push_back(a + shift);}
    }

    io::File file("temp/rigidbody/ausaxs_convsym_single_bad.pdb");
    write_pdb(file, assembly);

    Sequencer seq(io::ExistingFile("tests/files/SASDJG5.dat"));
    seq.setup().load(std::vector<std::string>{file});
    CHECK_THROWS([&]{ ConvertToSymmetryElement convert(&seq, {0}, "c3"); }());
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
