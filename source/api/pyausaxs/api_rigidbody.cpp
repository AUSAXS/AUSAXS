// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/pyausaxs/api_rigidbody.h>
#include <api/ObjectStorage.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <rigidbody/sequencer/detail/ValidElements.h>
#include <rigidbody/constraints/ConstrainedFitter.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>

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