// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/continuation/ContinuationState.h>
#include <rigidbody/constraints/AttractorConstraint.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/constraints/DistanceConstraintAtom.h>
#include <rigidbody/constraints/DistanceConstraintBond.h>
#include <rigidbody/constraints/DistanceConstraintCM.h>
#include <rigidbody/constraints/RepellerConstraint.h>
#include <rigidbody/Rigidbody.h>
#include <data/Body.h>
#include <data/Molecule.h>
#include <data/symmetry/CompositeSymmetry.h>
#include <data/symmetry/PredefinedSymmetries.h>
#include <data/symmetry/ReferenceSymmetry.h>
#include <io/File.h>
#include <utility/Exceptions.h>

#include <cstring>
#include <fstream>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::continuation;

namespace {
    // Format identity. The magic is checked so that pointing `load {continue ...}` at an unrelated file fails with a
    // clear message rather than a garbage read; the version is checked so that a future layout change is rejected
    // outright instead of being silently misread.
    constexpr char magic[16] = "AUSAXS-CONTINUE";
    constexpr std::uint32_t format_version = 1;

    // Tags distinguishing the three shapes an ISymmetry can take. A symmetry is written as a tree rather than a flat
    // record because a CompositeSymmetry owns two sub-symmetries, each with their own parameter spans.
    enum class SymmetryTag : std::uint8_t {leaf, composite, reference, reference_view};

    // Coordinates are always written as double regardless of what constants::coords_precision_t happens to be in this
    // build, so a state file stays readable across build configurations. Widening never loses anything.
    struct Writer {
        std::ofstream out;

        explicit Writer(const std::string& path) : out(path, std::ios::binary) {
            if (!out.is_open()) {throw except::io_error("ContinuationState: could not open \"" + path + "\" for writing.");}
        }

        template<typename T>
        void pod(const T& value) {out.write(reinterpret_cast<const char*>(&value), sizeof(T));}

        void u8(std::uint8_t v)   {pod(v);}
        void u32(std::uint32_t v) {pod(v);}
        void u64(std::uint64_t v) {pod(v);}
        void i32(std::int32_t v)  {pod(v);}
        void f64(double v)        {pod(v);}

        void str(const std::string& s) {
            u32(static_cast<std::uint32_t>(s.size()));
            out.write(s.data(), static_cast<std::streamsize>(s.size()));
        }

        void vec3(const Vector3<double>& v) {f64(v.x()); f64(v.y()); f64(v.z());}
    };

    struct Reader {
        std::ifstream in;
        std::string path;

        explicit Reader(const std::string& path) : in(path, std::ios::binary), path(path) {
            if (!in.is_open()) {throw except::io_error("ContinuationState: could not open \"" + path + "\" for reading.");}
        }

        template<typename T>
        T pod() {
            T value;
            in.read(reinterpret_cast<char*>(&value), sizeof(T));
            if (!in) {throw except::io_error("ContinuationState: \"" + path + "\" ended unexpectedly; the file is truncated or corrupt.");}
            return value;
        }

        std::uint8_t u8()   {return pod<std::uint8_t>();}
        std::uint32_t u32() {return pod<std::uint32_t>();}
        std::uint64_t u64() {return pod<std::uint64_t>();}
        std::int32_t i32()  {return pod<std::int32_t>();}
        double f64()        {return pod<double>();}

        std::string str() {
            auto size = u32();
            std::string s(size, '\0');
            in.read(s.data(), static_cast<std::streamsize>(size));
            if (!in) {throw except::io_error("ContinuationState: \"" + path + "\" ended unexpectedly; the file is truncated or corrupt.");}
            return s;
        }

        Vector3<double> vec3() {double x = f64(), y = f64(), z = f64(); return {x, y, z};}
    };

    // ----- symmetries ---------------------------------------------------------
    // A leaf is reconstructed by name (symmetry::create) and then has its optimisable parameters overwritten, which is
    // what keeps a refined pose intact: the name only fixes the *kind* of symmetry, the spans carry where it currently is.
    void write_symmetry(Writer& w, symmetry::ISymmetry& sym) {
        if (auto* composite = dynamic_cast<symmetry::CompositeSymmetry*>(&sym)) {
            w.u8(static_cast<std::uint8_t>(SymmetryTag::composite));
            write_symmetry(w, *composite->inner);
            write_symmetry(w, *composite->outer);
            return;
        }
        if (auto* view = dynamic_cast<symmetry::ReferenceSymmetryView*>(&sym)) {
            w.u8(static_cast<std::uint8_t>(SymmetryTag::reference_view));
            w.i32(view->primary_body);
            w.i32(view->symmetry_index);
            return;
        }
        if (auto* ref = dynamic_cast<symmetry::ReferenceSymmetry*>(&sym)) {
            w.u8(static_cast<std::uint8_t>(SymmetryTag::reference));
            w.u32(static_cast<std::uint32_t>(ref->bodies.size()));
            for (int body : ref->bodies) {w.i32(body);}
            w.u32(static_cast<std::uint32_t>(ref->slots.size()));
            for (int slot : ref->slots) {w.i32(slot);}
            write_symmetry(w, *ref->base);
            return;
        }
        w.u8(static_cast<std::uint8_t>(SymmetryTag::leaf));
        w.str(sym.type_name());
        auto translation = sym.span_translation();
        w.u32(static_cast<std::uint32_t>(translation.size()));
        for (double v : translation) {w.f64(v);}
        auto rotation = sym.span_rotation();
        w.u32(static_cast<std::uint32_t>(rotation.size()));
        for (double v : rotation) {w.f64(v);}
    }

    std::unique_ptr<symmetry::ISymmetry> read_symmetry(Reader& r, observer_ptr<const data::Molecule> molecule) {
        switch (static_cast<SymmetryTag>(r.u8())) {
            case SymmetryTag::composite: {
                auto inner = read_symmetry(r, molecule);
                auto outer = read_symmetry(r, molecule);
                return std::make_unique<symmetry::CompositeSymmetry>(std::move(inner), std::move(outer));
            }
            case SymmetryTag::reference_view: {
                int primary_body = r.i32();
                int symmetry_index = r.i32();
                return std::make_unique<symmetry::ReferenceSymmetryView>(molecule, primary_body, symmetry_index);
            }
            case SymmetryTag::reference: {
                std::vector<int> bodies(r.u32());
                for (auto& body : bodies) {body = r.i32();}
                std::vector<int> slots(r.u32());
                for (auto& slot : slots) {slot = r.i32();}
                auto base = read_symmetry(r, molecule);
                return std::make_unique<symmetry::ReferenceSymmetry>(std::move(base), std::move(bodies), std::move(slots), molecule);
            }
            case SymmetryTag::leaf: {
                auto name = r.str();
                auto sym = symmetry::create(name);
                auto translation = sym->span_translation();
                auto n_translation = r.u32();
                if (n_translation != translation.size()) {
                    throw except::io_error(
                        "ContinuationState: symmetry \"" + name + "\" expects " + std::to_string(translation.size()) +
                        " translation parameters, but the file holds " + std::to_string(n_translation) + "."
                    );
                }
                for (auto& v : translation) {v = r.f64();}
                auto rotation = sym->span_rotation();
                auto n_rotation = r.u32();
                if (n_rotation != rotation.size()) {
                    throw except::io_error(
                        "ContinuationState: symmetry \"" + name + "\" expects " + std::to_string(rotation.size()) +
                        " rotation parameters, but the file holds " + std::to_string(n_rotation) + "."
                    );
                }
                for (auto& v : rotation) {v = r.f64();}
                return sym;
            }
        }
        throw except::io_error("ContinuationState: unknown symmetry tag; the file is corrupt.");
    }

    // Advance past one symmetry record without building anything. Used by the first read pass, which runs before the
    // molecule exists and so cannot construct the reference symmetries that need to point at it.
    void skip_symmetry(Reader& r) {
        switch (static_cast<SymmetryTag>(r.u8())) {
            case SymmetryTag::composite:
                skip_symmetry(r);
                skip_symmetry(r);
                return;
            case SymmetryTag::reference_view:
                r.i32(); r.i32();
                return;
            case SymmetryTag::reference: {
                for (auto n = r.u32(); n > 0; --n) {r.i32();}
                for (auto n = r.u32(); n > 0; --n) {r.i32();}
                skip_symmetry(r);
                return;
            }
            case SymmetryTag::leaf: {
                r.str();
                for (auto n = r.u32(); n > 0; --n) {r.f64();}
                for (auto n = r.u32(); n > 0; --n) {r.f64();}
                return;
            }
        }
        throw except::io_error("ContinuationState: unknown symmetry tag; the file is corrupt.");
    }

    // ----- constraints --------------------------------------------------------
    // The most-derived type is identified first, since every one of these is a base of the next. An unrecognised
    // subclass throws rather than falling back to `atom`: silently writing a constraint out as a weaker kind would
    // change what the continued run is actually optimising.
    ContinuationConstraint::Kind classify(const constraints::IDistanceConstraint& c) {
        if (dynamic_cast<const constraints::AttractorConstraint*>(&c))    {return ContinuationConstraint::Kind::attractor;}
        if (dynamic_cast<const constraints::RepellerConstraint*>(&c))     {return ContinuationConstraint::Kind::repeller;}
        if (dynamic_cast<const constraints::DistanceConstraintCM*>(&c))   {return ContinuationConstraint::Kind::cm;}
        if (dynamic_cast<const constraints::DistanceConstraintBond*>(&c)) {return ContinuationConstraint::Kind::bond;}
        if (dynamic_cast<const constraints::DistanceConstraintAtom*>(&c)) {return ContinuationConstraint::Kind::atom;}
        throw except::io_error("ContinuationState: cannot store an unrecognised distance constraint type.");
    }

    // Every distance constraint the manager holds, from both of its lists. The split between them is about how
    // constraints are *looked up* during refinement (only backbone bonds are indexed per body), not about what they
    // are, so a writer that read only the discoverable list would quietly drop every cm/attract/repel constraint.
    // Constraints that are not distance constraints - notably the overlap constraint, which every Rigidbody creates
    // for itself - are left out and simply regenerated by the continued run.
    std::vector<observer_ptr<const constraints::IDistanceConstraint>> distance_constraints(const constraints::ConstraintManager& manager) {
        std::vector<observer_ptr<const constraints::IDistanceConstraint>> out;
        for (const auto& c : manager.discoverable_constraints) {out.push_back(c.get());}
        for (const auto& c : manager.non_discoverable_constraints) {
            if (auto* distance = dynamic_cast<const constraints::IDistanceConstraint*>(c.get())) {out.push_back(distance);}
        }
        return out;
    }
}

void ausaxs::rigidbody::continuation::write_continuation_state(
    const ausaxs::io::File& path, const Rigidbody& rigidbody, const std::vector<std::string>& body_names
) {
    path.directory().create();
    Writer w(path.str());

    w.out.write(magic, sizeof(magic));
    w.u32(format_version);

    const auto& bodies = rigidbody.molecule.get_bodies();
    w.u32(static_cast<std::uint32_t>(bodies.size()));
    w.u32(static_cast<std::uint32_t>(body_names.size()));
    for (const auto& name : body_names) {w.str(name);}

    for (const auto& body : bodies) {
        const auto& atoms = body.get_atoms();
        w.u64(atoms.size());
        for (const auto& atom : atoms) {
            w.vec3(atom.coordinates());
            w.f64(atom.weight());
            w.i32(static_cast<std::int32_t>(atom.form_factor_type()));
        }

        // per-atom metadata is what a .pdb round-trip destroys: without the residue ids and Cα flags, a continued run
        // could neither re-split by residue nor generate backbone constraints.
        const auto& metadata = body.get_metadata();
        std::uint8_t flags = 0;
        if (metadata) {
            if (metadata->backbone)    {flags |= 1 << 0;}
            if (metadata->residue_seq) {flags |= 1 << 1;}
            if (metadata->occupancy)   {flags |= 1 << 2;}
        }
        w.u8(flags);
        if (flags & (1 << 0)) {for (auto v : *metadata->backbone)    {w.u8(static_cast<std::uint8_t>(v));}}
        if (flags & (1 << 1)) {for (auto v : *metadata->residue_seq) {w.i32(v);}}
        if (flags & (1 << 2)) {for (auto v : *metadata->occupancy)   {w.f64(v);}}

        // Hydration is deliberately not stored. It is derived data, not part of the rigid-body state: every refinement
        // calls generate_new_hydration(), which clears the layer and rebuilds it from the grid. Restoring the previous
        // run's shell would not survive that, but it *would* be present while the grid is first built, perturbing the
        // shell that gets regenerated and so shifting chi2 for a structure that has not moved at all.

        // read through the const facade so the bodies aren't flagged as modified merely by being saved; the const_cast is
        // only needed because span_translation()/span_rotation() are non-const accessors on an otherwise read-only walk.
        const auto& symmetries = body.symmetry().get();
        w.u32(static_cast<std::uint32_t>(symmetries.size()));
        for (const auto& sym : symmetries) {write_symmetry(w, const_cast<symmetry::ISymmetry&>(*sym));}
    }

    auto constraints = distance_constraints(*rigidbody.constraints);
    w.u32(static_cast<std::uint32_t>(constraints.size()));
    for (const auto& c : constraints) {
        w.u8(static_cast<std::uint8_t>(classify(*c)));
        w.i32(c->ibody1); w.i32(c->iatom1);
        w.i32(c->ibody2); w.i32(c->iatom2);
        w.i32(c->isym1.first); w.i32(c->isym1.second);
        w.i32(c->isym2.first); w.i32(c->isym2.second);
        w.f64(c->d_target);
    }
}

ContinuationState ausaxs::rigidbody::continuation::read_continuation_state(const ausaxs::io::File& path) {
    Reader r(path.str());

    char header[sizeof(magic)];
    r.in.read(header, sizeof(header));
    if (!r.in || std::memcmp(header, magic, sizeof(magic)) != 0) {
        throw except::io_error("ContinuationState: \"" + path.str() + "\" is not a continuation state file.");
    }
    if (auto version = r.u32(); version != format_version) {
        throw except::io_error(
            "ContinuationState: \"" + path.str() + "\" was written by format version " + std::to_string(version) +
            ", but this build reads version " + std::to_string(format_version) + "."
        );
    }

    auto n_bodies = r.u32();
    ContinuationState state;
    state.body_names.resize(r.u32());
    for (auto& name : state.body_names) {name = r.str();}

    // Bodies are built without their symmetries first: a ReferenceSymmetry holds a pointer to the molecule, so the
    // molecule has to exist before any symmetry can be constructed. The raw symmetry records are parked until then.
    std::vector<data::Body> bodies;
    bodies.reserve(n_bodies);
    std::vector<std::streampos> symmetry_offsets;
    symmetry_offsets.reserve(n_bodies);

    for (std::uint32_t i = 0; i < n_bodies; ++i) {
        auto n_atoms = r.u64();
        std::vector<data::AtomFF> atoms;
        atoms.reserve(n_atoms);
        for (std::uint64_t j = 0; j < n_atoms; ++j) {
            auto coords = r.vec3();
            auto weight = r.f64();
            auto type = static_cast<form_factor::form_factor_t>(r.i32());
            atoms.emplace_back(coords, type, weight);
        }

        data::AtomMetadata metadata;
        auto flags = r.u8();
        if (flags & (1 << 0)) {
            std::vector<data::backbone_t> backbone(n_atoms);
            for (auto& v : backbone) {v = static_cast<data::backbone_t>(r.u8());}
            metadata.backbone = std::move(backbone);
        }
        if (flags & (1 << 1)) {
            std::vector<int> residue_seq(n_atoms);
            for (auto& v : residue_seq) {v = r.i32();}
            metadata.residue_seq = std::move(residue_seq);
        }
        if (flags & (1 << 2)) {
            std::vector<float> occupancy(n_atoms);
            for (auto& v : occupancy) {v = static_cast<float>(r.f64());}
            metadata.occupancy = std::move(occupancy);
        }

        data::Body body(std::move(atoms));
        if (flags) {body.set_metadata(std::move(metadata));}
        bodies.push_back(std::move(body));

        // skip past this body's symmetries; they are read in a second pass once the molecule exists
        symmetry_offsets.push_back(r.in.tellg());
        auto n_symmetries = r.u32();
        for (std::uint32_t j = 0; j < n_symmetries; ++j) {skip_symmetry(r);}
    }

    auto constraints_offset = r.in.tellg();
    state.molecule = data::Molecule(std::move(bodies));

    for (std::uint32_t i = 0; i < n_bodies; ++i) {
        r.in.seekg(symmetry_offsets[i]);
        auto n_symmetries = r.u32();
        for (std::uint32_t j = 0; j < n_symmetries; ++j) {
            state.molecule.get_body(i).symmetry().add(read_symmetry(r, &state.molecule));
        }
    }

    r.in.seekg(constraints_offset);
    state.constraints.resize(r.u32());
    for (auto& c : state.constraints) {
        c.kind = static_cast<ContinuationConstraint::Kind>(r.u8());
        c.ibody1 = r.i32(); c.iatom1 = r.i32();
        c.ibody2 = r.i32(); c.iatom2 = r.i32();
        c.isym1.first = r.i32(); c.isym1.second = r.i32();
        c.isym2.first = r.i32(); c.isym2.second = r.i32();
        c.d_target = r.f64();
    }
    return state;
}
