// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/pyausaxs/api_rigidbody.h>
#include <api/ObjectStorage.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/elements/setup/LoadElementWithMetadata.h>
#include <rigidbody/sequencer/elements/UpdateElement.h>
#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/detail/ValidElements.h>
#include <rigidbody/constraints/ConstrainedFitter.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <data/Molecule.h>
#include <data/Body.h>
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

namespace {
    // Swap the load element keyword (`load`/`open`, as the first token on a line) for the hidden
    // `load_preview` keyword, so the structure is loaded with the per-atom metadata previews need.
    std::string substitute_load_with_metadata(const std::string& script) {
        std::string result;
        result.reserve(script.size() + 32);
        std::size_t i = 0;
        while (i < script.size()) {
            std::size_t eol = script.find('\n', i);
            std::size_t line_end = (eol == std::string::npos) ? script.size() : eol;

            std::size_t token_begin = i;
            while (token_begin < line_end && (script[token_begin] == ' ' || script[token_begin] == '\t')) {++token_begin;}
            std::size_t token_end = token_begin;
            while (token_end < line_end && script[token_end] != ' ' && script[token_end] != '\t' && script[token_end] != '{') {++token_end;}

            result.append(script, i, token_begin - i); // leading whitespace
            std::string token = script.substr(token_begin, token_end - token_begin);
            result += (token == "load" || token == "open") ? "load_preview" : token;
            result.append(script, token_end, line_end - token_end); // remainder of the line
            if (eol != std::string::npos) {result += '\n';}
            i = line_end + 1;
            if (eol == std::string::npos) {break;}
        }
        return result;
    }
}

struct _rigidbody_preview_structure_obj {
    std::vector<double> x, y, z;
    std::vector<int> body_index, copy_index, residue_seq, is_ca;
};
int rigidbody_get_preview_structure(
    int rigidbody_id,
    double** x, double** y, double** z,
    int** body_index, int** copy_index, int** residue_seq, int** is_ca,
    int* n_atoms,
    int* status
) {return execute_with_catch([&]() {
    auto script_obj = api::ObjectStorage::get_object<_rigidbody_script_obj>(rigidbody_id);
    if (!script_obj) {ErrorMessage::last_error = "Invalid rigidbody script id: \"" + std::to_string(rigidbody_id) + "\""; return -1;}

    auto sequencer = rigidbody::sequencer::SequenceParser().parse_text(substitute_load_with_metadata(script_obj->script));
    auto molecule = sequencer->_get_molecule();
    // the load_preview element refreshes these statics during the parse above; aligned with the
    // molecule's base atoms in body order (i.e. before symmetry copies are realized)
    const auto& meta_res = rigidbody::sequencer::LoadElementWithMetadata::residue_seq;
    const auto& meta_ca  = rigidbody::sequencer::LoadElementWithMetadata::is_ca;
    bool has_meta = !meta_res.empty();

    _rigidbody_preview_structure_obj data;
    int base_offset = 0; // running index into the per-base-atom metadata
    int bidx = 0;
    for (const auto& body : molecule->get_bodies()) {
        int na = static_cast<int>(body.size_atom());

        // explicit structure of this body: na*(1 + n_copies) atoms, laid out as
        // [original, copy_1, copy_2, ...], each block reusing the base atom order
        auto bstruct = body.symmetry().explicit_structure();
        int blocks = (na > 0) ? static_cast<int>(bstruct.atoms.size()) / na : 0;
        for (int copy = 0; copy < blocks; ++copy) {
            for (int i = 0; i < na; ++i) {
                const auto& atom = bstruct.atoms[copy*na + i];
                data.x.push_back(atom.x());
                data.y.push_back(atom.y());
                data.z.push_back(atom.z());
                data.body_index.push_back(bidx);
                data.copy_index.push_back(copy);
                int gi = base_offset + i; // copies share their original atom's metadata
                data.residue_seq.push_back(has_meta && gi < static_cast<int>(meta_res.size()) ? meta_res[gi] : -1);
                data.is_ca.push_back(has_meta && gi < static_cast<int>(meta_ca.size()) ? static_cast<int>(meta_ca[gi]) : 0);
            }
        }
        base_offset += na;
        ++bidx;
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