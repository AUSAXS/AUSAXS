// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/pyausaxs/api_rigidbody.h>
#include <api/ObjectStorage.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/elements/UpdateElement.h>
#include <rigidbody/sequencer/elements/LoopElement.h>
#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/detail/ValidElements.h>
#include <rigidbody/constraints/ConstrainedFitter.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/constraints/AttractorConstraint.h>
#include <rigidbody/constraints/RepellerConstraint.h>
#include <rigidbody/constraints/DistanceConstraintCM.h>
#include <rigidbody/constraints/IDistanceConstraint.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <data/atoms/AtomMetadata.h>
#include <data/symmetry/BodySymmetryFacade.h>
#include <data/detail/SimpleBody.h>
#include <settings/MoleculeSettings.h>
#include <utility/Exceptions.h>

using namespace ausaxs;

namespace {
struct _rigidbody_script_obj {
    std::string script;
    // lazily parsed and cached on the first read-only preview query (see get_cached_sequencer below)
    std::unique_ptr<rigidbody::sequencer::Sequencer> sequencer;
};
}
int rigidbody_load_script(
    const char* script,
    int* status
) {return execute_with_catch([&]() {
    _rigidbody_script_obj data;
    data.script = script;
    int data_id = api::ObjectStorage::register_object(std::move(data));
    return data_id;
}, status);}

namespace {
    rigidbody::sequencer::Sequencer& get_cached_sequencer(_rigidbody_script_obj& obj) {
        if (!obj.sequencer) {
            // every current caller of the cached path needs Cα/backbone metadata for its output
            settings::molecule::store_calpha = true;
            obj.sequencer = rigidbody::sequencer::SequenceParser().parse_text(obj.script);
        }
        return *obj.sequencer;
    }
}

namespace {
struct _rigidbody_preview_structure_obj {
    std::vector<double> x, y, z;
    std::vector<int> body_index, copy_index, residue_seq, is_ca;
    std::vector<int> constraint_data; // flat triplets: [index1, index2, type, ...]
};
}
int rigidbody_get_preview_structure(
    int rigidbody_id,
    double** x, double** y, double** z,
    int** body_index, int** copy_index, int** residue_seq, int** is_ca,
    int* n_atoms,
    int** constraint_data, int* n_constraints,
    int* status
) {return execute_with_catch([&]() {
    auto script_obj = api::ObjectStorage::get_object<_rigidbody_script_obj>(rigidbody_id);
    if (!script_obj) {ErrorMessage::last_error = "Invalid rigidbody script id: \"" + std::to_string(rigidbody_id) + "\""; return -1;}

    auto& sequencer = get_cached_sequencer(*script_obj);
    auto molecule = sequencer._get_molecule();

    _rigidbody_preview_structure_obj data;
    int bidx = 0;
    // flat index of atom 0 of each body's copy 0; used below to map constraint atom indices
    std::vector<int> body_atom0_starts;
    for (const auto& body : molecule->get_bodies()) {
        int na = static_cast<int>(body.size_atom());

        // each body's metadata is parallel-indexed to its own base atoms (symmetry copies reuse it)
        const auto& md = body.get_metadata();
        const std::vector<int>* res_seq = (md && md->residue_seq) ? &md->residue_seq.value() : nullptr;
        const std::vector<data::backbone_t>* backbone = (md && md->backbone) ? &md->backbone.value() : nullptr;

        // explicit structure of this body: na*(1 + n_copies) atoms, laid out as
        // [original, copy_1, copy_2, ...], each block reusing the base atom order
        auto bstruct = body.symmetry().explicit_structure();
        int blocks = (na > 0) ? static_cast<int>(bstruct.atoms.size()) / na : 0;
        body_atom0_starts.push_back(static_cast<int>(data.x.size())); // copy 0 starts here
        for (int copy = 0; copy < blocks; ++copy) {
            for (int i = 0; i < na; ++i) {
                const auto& atom = bstruct.atoms[copy*na + i];
                data.x.push_back(atom.x());
                data.y.push_back(atom.y());
                data.z.push_back(atom.z());
                data.body_index.push_back(bidx);
                data.copy_index.push_back(copy);
                data.residue_seq.push_back(res_seq ? (*res_seq)[i] : -1);
                data.is_ca.push_back(backbone && (*backbone)[i] == data::backbone_t::c_alpha ? 1 : 0);
            }
        }
        ++bidx;
    }

    // constraint type codes: 0=backbone, 1=CM, 2=attractor, 3=repulsor
    // user-generated constraints always use the base body (isym={-1,-1}, copy=0), so iatom directly indexes into copy 0's atom block
    auto emit_constraint = [&] (const rigidbody::constraints::IDistanceConstraint* c) {
        int idx1 = body_atom0_starts[c->ibody1] + c->iatom1;
        int idx2 = body_atom0_starts[c->ibody2] + c->iatom2;
        int type;
        if      (dynamic_cast<const rigidbody::constraints::AttractorConstraint*>(c)) {type = 2;}
        else if (dynamic_cast<const rigidbody::constraints::RepellerConstraint*> (c)) {type = 3;}
        else if (dynamic_cast<const rigidbody::constraints::DistanceConstraintCM*>(c)) {type = 1;}
        else {type = 0;}
        data.constraint_data.push_back(idx1);
        data.constraint_data.push_back(idx2);
        data.constraint_data.push_back(type);
    };
    for (const auto& c : sequencer._get_rigidbody()->constraints->discoverable_constraints) {
        emit_constraint(c.get());
    }
    for (const auto& c : sequencer._get_rigidbody()->constraints->non_discoverable_constraints) {
        if (auto* dc = dynamic_cast<const rigidbody::constraints::IDistanceConstraint*>(c.get())) {
            emit_constraint(dc);
        }
    }

    int data_id = api::ObjectStorage::register_object(std::move(data));
    auto ref = api::ObjectStorage::get_object<_rigidbody_preview_structure_obj>(data_id);
    *x = ref->x.data();
    *y = ref->y.data();
    *z = ref->z.data();
    *body_index = ref->body_index.data();
    *copy_index = ref->copy_index.data();
    *residue_seq = ref->residue_seq.data();
    *is_ca = ref->is_ca.data();
    *n_atoms = static_cast<int>(ref->x.size());
    *constraint_data = ref->constraint_data.empty() ? nullptr : ref->constraint_data.data();
    *n_constraints = static_cast<int>(ref->constraint_data.size()) / 3;
    return data_id;
}, status);}

namespace {
struct _rigidbody_live_structure_obj {
    std::vector<double> x, y, z;
};
}
int rigidbody_get_live_structure(
    double** x, double** y, double** z,
    int* n_atoms, int* version,
    int* status
) {return execute_with_catch([&]() {
    _rigidbody_live_structure_obj data;
    int ver = 0;

    rigidbody::sequencer::UpdateElement::lock();
    data.x = rigidbody::sequencer::UpdateElement::x;
    data.y = rigidbody::sequencer::UpdateElement::y;
    data.z = rigidbody::sequencer::UpdateElement::z;
    ver = rigidbody::sequencer::UpdateElement::version;
    rigidbody::sequencer::UpdateElement::unlock();

    int data_id = api::ObjectStorage::register_object(std::move(data));
    auto ref = api::ObjectStorage::get_object<_rigidbody_live_structure_obj>(data_id);
    *x = ref->x.data();
    *y = ref->y.data();
    *z = ref->z.data();
    *n_atoms = static_cast<int>(ref->x.size());
    *version = ver;
    return data_id;
}, status);}

namespace {
struct _rigidbody_live_poller_obj {
    std::vector<double> x, y, z;
};
}
int rigidbody_create_live_poller(
    int* status
) {return execute_with_catch([&]() {
    return api::ObjectStorage::register_object(_rigidbody_live_poller_obj{});
}, status);}

void rigidbody_poll_live_structure(
    int poller_id,
    double** x, double** y, double** z,
    int* n_atoms, int* version,
    int* status
) {execute_with_catch([&]() {
    auto poller = api::ObjectStorage::get_object<_rigidbody_live_poller_obj>(poller_id);
    if (!poller) {throw except::invalid_argument("Invalid live poller id: \"" + std::to_string(poller_id) + "\"");}

    // copied under the publisher's lock: the buffers are resized from a thread pool task, so reading them
    // directly from the consumer's thread would race with a reallocation
    rigidbody::sequencer::UpdateElement::lock();
    poller->x = rigidbody::sequencer::UpdateElement::x;
    poller->y = rigidbody::sequencer::UpdateElement::y;
    poller->z = rigidbody::sequencer::UpdateElement::z;
    *version = rigidbody::sequencer::UpdateElement::version;
    rigidbody::sequencer::UpdateElement::unlock();

    *x = poller->x.data();
    *y = poller->y.data();
    *z = poller->z.data();
    *n_atoms = static_cast<int>(poller->x.size());
}, status);}

void rigidbody_register_live_consumer(bool connected, int* status) {execute_with_catch([&]() {
    rigidbody::sequencer::UpdateElement::live_consumer_connected = connected;
}, status);}

void rigidbody_validate(
    int rigidbody_id,
    int* status
) {return execute_with_catch([&]() {
    auto script_obj = api::ObjectStorage::get_object<_rigidbody_script_obj>(rigidbody_id);
    if (!script_obj) {ErrorMessage::last_error = "Invalid rigidbody script id: \"" + std::to_string(rigidbody_id) + "\""; return;}
    rigidbody::sequencer::SequenceParser().parse_text(script_obj->script);
}, status);}

namespace {
struct _rigidbody_get_data_obj {
    std::vector<double> q, I, I_err, I_inter;
};
}
int rigidbody_run(
    int rigidbody_id,
    double** q, double** I, double** I_err, double** I_interp, int* n_points,
    int* status
) {return execute_with_catch([&]() {
    auto script_obj = api::ObjectStorage::get_object<_rigidbody_script_obj>(rigidbody_id);
    if (!script_obj) {ErrorMessage::last_error = "Invalid rigidbody script id: \"" + std::to_string(rigidbody_id) + "\""; return -1;}
    auto sequencer = rigidbody::sequencer::SequenceParser().parse_text(script_obj->script);
    sequencer->execute();

    auto fit_result = sequencer->_get_controller()->get_fitter()->fit();
    auto data = fit_result->curves.select_columns({0, 1, 2, 3});
    _rigidbody_get_data_obj data_obj;
    data_obj.q = data.col(0);
    data_obj.I = data.col(1);
    data_obj.I_err = data.col(2);
    data_obj.I_inter = data.col(3);
    int data_id = api::ObjectStorage::register_object(std::move(data_obj));
    auto ref = api::ObjectStorage::get_object<_rigidbody_get_data_obj>(data_id);
    *q = ref->q.data();
    *I = ref->I.data();
    *I_err = ref->I_err.data();
    *I_interp = ref->I_inter.data();
    *n_points = static_cast<int>(ref->q.size());
    return data_id;
}, status);}

void rigidbody_stop_run(
    int* status
) {execute_with_catch([&]() {
    rigidbody::sequencer::LoopElement::_request_stop();
}, status);}

void rigidbody_get_valid_elements(
    const char*** elements,
    int* size,
    int* status
) {execute_with_catch([&]() {
    static std::vector<std::string> valid_elements = rigidbody::sequencer::detail::valid_elements();
    static std::vector<const char*> valid_elements_cstr = [&] () {
        std::vector<const char*> cstrs;
        for (const auto& elem : valid_elements) {cstrs.push_back(elem.c_str());}
        return cstrs;
    }();
    *elements = valid_elements_cstr.data();
    *size = valid_elements.size();
}, status);}

void rigidbody_get_valid_arguments(
    const char* element_name,
    const char*** arguments,
    int* size,
    int* status
) {execute_with_catch([&]() {
    auto type = rigidbody::sequencer::detail::get_type(element_name);
    static std::unordered_map<rigidbody::sequencer::detail::ElementType, std::vector<std::string>> valid_arguments_map;
    static std::unordered_map<rigidbody::sequencer::detail::ElementType, std::vector<const char*>> valid_arguments_cstr_map;
    if (!valid_arguments_map.contains(type)) {
        valid_arguments_map[type] = rigidbody::sequencer::detail::valid_arguments(type);
        std::vector<const char*> cstrs;
        for (const auto& arg : valid_arguments_map[type]) {cstrs.push_back(arg.c_str());}
        valid_arguments_cstr_map[type] = cstrs;
    }
    *arguments = valid_arguments_cstr_map[type].data();
    *size = valid_arguments_map[type].size();
}, status);}

void rigidbody_get_valid_inline_arguments(
    const char* element_name,
    const char*** arguments,
    int* size,
    int* min_count,
    int* max_count,
    int* status
) {execute_with_catch([&]() {
    auto type = rigidbody::sequencer::detail::get_type(element_name);
    static std::unordered_map<rigidbody::sequencer::detail::ElementType, rigidbody::sequencer::InlineSignature> signature_map;
    static std::unordered_map<rigidbody::sequencer::detail::ElementType, std::vector<const char*>> signature_cstr_map;
    if (!signature_map.contains(type)) {
        signature_map[type] = rigidbody::sequencer::detail::valid_inline_arguments(type);
        std::vector<const char*> cstrs;
        for (const auto& arg : signature_map[type].names) {cstrs.push_back(arg.c_str());}
        signature_cstr_map[type] = cstrs;
    }
    *arguments = signature_cstr_map[type].data();
    *size = signature_map[type].names.size();
    *min_count = static_cast<int>(signature_map[type].min);
    *max_count = static_cast<int>(signature_map[type].max);
}, status);}

void rigidbody_get_body_names(
    int rigidbody_id,
    const char*** names,
    int* size,
    int* status
) {execute_with_catch([&]() {
    auto script_obj = api::ObjectStorage::get_object<_rigidbody_script_obj>(rigidbody_id);
    if (!script_obj) {ErrorMessage::last_error = "Invalid rigidbody script id: \"" + std::to_string(rigidbody_id) + "\""; return;}

    auto& sequencer = get_cached_sequencer(*script_obj);
    // the setup elements (merge/delete/convert_to_symmetry) are applied during parsing, so the registry
    // already reflects the final body set, ordered identically to rigidbody_get_preview_structure's bodies
    static std::vector<std::string> body_names;
    static std::vector<const char*> body_names_cstr;
    body_names = sequencer.setup()._body_name_registry().base_body_names();
    body_names_cstr.clear();
    body_names_cstr.reserve(body_names.size());
    for (const auto& name : body_names) {body_names_cstr.push_back(name.c_str());}

    *names = body_names_cstr.data();
    *size = static_cast<int>(body_names_cstr.size());
}, status);}

namespace {
struct _rigidbody_symmetry_layout_obj {
    std::vector<int> body, copy, symmetry, replica;
    std::vector<std::string> type, name;
    std::vector<const char*> type_ptr, name_ptr;
};
}
int rigidbody_get_symmetry_layout(
    int rigidbody_id,
    int** body, int** copy, int** symmetry, int** replica,
    const char*** type, const char*** name,
    int* n_replicas,
    int* status
) {return execute_with_catch([&]() {
    auto script_obj = api::ObjectStorage::get_object<_rigidbody_script_obj>(rigidbody_id);
    if (!script_obj) {ErrorMessage::last_error = "Invalid rigidbody script id: \"" + std::to_string(rigidbody_id) + "\""; return -1;}

    auto& sequencer = get_cached_sequencer(*script_obj);
    auto molecule = sequencer._get_molecule();
    const auto& name_registry = sequencer.setup()._body_name_registry();

    _rigidbody_symmetry_layout_obj data;
    int bidx = 0;
    for (const auto& body_obj : molecule->get_bodies()) {
        // copy 0 is the unmodified original; the remaining copies are laid out sequentially per
        // symmetry, exactly as explicit_structure() (and thus rigidbody_get_preview_structure) builds them
        int copy_idx = 1;
        int isymmetry = 0;
        for (const auto& sym_ptr : body_obj.symmetry().get()) {
            int reps = static_cast<int>(sym_ptr->repetitions());
            std::string type_name = sym_ptr->type_name();
            for (int replica_idx = 1; replica_idx <= reps; ++replica_idx) {
                data.body.push_back(bidx);
                data.copy.push_back(copy_idx);
                data.symmetry.push_back(isymmetry);
                data.replica.push_back(replica_idx);
                data.type.push_back(type_name);
                data.name.push_back(name_registry.entry(rigidbody::sequencer::detail::to_index(bidx, isymmetry, replica_idx)).display_name());
                ++copy_idx;
            }
            ++isymmetry;
        }
        ++bidx;
    }

    data.type_ptr.reserve(data.type.size());
    for (const auto& s : data.type) {data.type_ptr.push_back(s.c_str());}
    data.name_ptr.reserve(data.name.size());
    for (const auto& s : data.name) {data.name_ptr.push_back(s.c_str());}

    int data_id = api::ObjectStorage::register_object(std::move(data));
    auto ref = api::ObjectStorage::get_object<_rigidbody_symmetry_layout_obj>(data_id);
    *body = ref->body.empty() ? nullptr : ref->body.data();
    *copy = ref->copy.empty() ? nullptr : ref->copy.data();
    *symmetry = ref->symmetry.empty() ? nullptr : ref->symmetry.data();
    *replica = ref->replica.empty() ? nullptr : ref->replica.data();
    *type = ref->type_ptr.empty() ? nullptr : ref->type_ptr.data();
    *name = ref->name_ptr.empty() ? nullptr : ref->name_ptr.data();
    *n_replicas = static_cast<int>(ref->body.size());
    return data_id;
}, status);}