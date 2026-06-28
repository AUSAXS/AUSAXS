// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/pyausaxs/api_rigidbody.h>
#include <api/ObjectStorage.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/elements/UpdateElement.h>
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
#include <settings/MoleculeSettings.h>
#include <data/symmetry/BodySymmetryFacade.h>
#include <data/detail/SimpleBody.h>

using namespace ausaxs;

struct _rigidbody_script_obj {
    std::string script;
};
int rigidbody_load_script(
    const char* script,
    int* status
) {return execute_with_catch([&]() {
    _rigidbody_script_obj data;
    data.script = script;
    int data_id = api::ObjectStorage::register_object(std::move(data));
    return data_id;
}, status);}

struct _rigidbody_preview_structure_obj {
    std::vector<double> x, y, z;
    std::vector<int> body_index, copy_index, residue_seq, is_ca;
    std::vector<int> constraint_data; // flat triplets: [index1, index2, type, ...]
};
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

    // enable per-atom metadata capture so the loaded bodies carry the residue id + Cα flag the
    // preview needs; restored afterwards to avoid leaking the setting into later API calls
    bool prev_store_calpha = settings::molecule::store_calpha;
    settings::molecule::store_calpha = true;
    std::unique_ptr<rigidbody::sequencer::Sequencer> sequencer;
    try {
        sequencer = rigidbody::sequencer::SequenceParser().parse_text(script_obj->script);
    } catch (...) {
        settings::molecule::store_calpha = prev_store_calpha;
        throw;
    }
    settings::molecule::store_calpha = prev_store_calpha;
    auto molecule = sequencer->_get_molecule();

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
    for (const auto& c : sequencer->_get_rigidbody()->constraints->discoverable_constraints) {
        emit_constraint(c.get());
    }
    for (const auto& c : sequencer->_get_rigidbody()->constraints->non_discoverable_constraints) {
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

struct _rigidbody_live_structure_obj {
    std::vector<double> x, y, z;
};
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

struct _data_get_data_obj {
    std::vector<double> q, I, I_err, I_inter;
};
int rigidbody_run(
    int rigidbody_id,
    double** q, double** I, double** I_err, double** I_interp, int* n_points,
    int* status
) {return execute_with_catch([&]() {
    auto script_obj = api::ObjectStorage::get_object<_rigidbody_script_obj>(rigidbody_id);
    if (!script_obj) {ErrorMessage::last_error = "Invalid rigidbody script id: \"" + std::to_string(rigidbody_id) + "\""; return -1;}
    auto sequencer = rigidbody::sequencer::SequenceParser().parse_text(script_obj->script);
    sequencer->execute();

    auto data = sequencer->_get_controller()->get_fitter()->fit()->curves.select_columns({0, 1, 2, 3});
    _data_get_data_obj data_obj;
    data_obj.q = data.col(0);
    data_obj.I = data.col(1);
    data_obj.I_err = data.col(2);
    data_obj.I_inter = data.col(3);
    int data_id = api::ObjectStorage::register_object(std::move(data_obj));
    auto ref = api::ObjectStorage::get_object<_data_get_data_obj>(data_id);
    *q = ref->q.data();
    *I = ref->I.data();
    *I_err = ref->I_err.data();
    *I_interp = ref->I_inter.data();
    *n_points = static_cast<int>(ref->q.size());
    return data_id;
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