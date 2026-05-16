// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/pyausaxs/api_pdb.h>
#include <api/ObjectStorage.h>
#include <rigidbody/sequencer/detail/SequenceParser.h>

using namespace ausaxs;

struct _rigidbody_script_obj {
    std::string script;
};
extern "C" API int rigidbody_load_script(
    const char* script,
    int* status
) {return execute_with_catch([&]() {
    _rigidbody_script_obj data;
    data.script = script;
    int data_id = api::ObjectStorage::register_object(std::move(data));
    return data_id;
}, status);}

extern "C" API int rigidbody_validate(
    int rigidbody_id,
    int* status
) {return execute_with_catch([&]() {
    auto script_obj = api::ObjectStorage::get_object<_rigidbody_script_obj>(rigidbody_id);
    if (!script_obj) {ErrorMessage::last_error = "Invalid rigidbody script id: \"" + std::to_string(rigidbody_id) + "\""; return -1;}
    rigidbody::sequencer::SequenceParser().parse_text(script_obj->script);
    return 0;
}, status);}

extern "C" API int rigidbody_run(
    int rigidbody_id,
    int* status
) {return execute_with_catch([&]() {
    auto script_obj = api::ObjectStorage::get_object<_rigidbody_script_obj>(rigidbody_id);
    if (!script_obj) {ErrorMessage::last_error = "Invalid rigidbody script id: \"" + std::to_string(rigidbody_id) + "\""; return -1;}
    auto sequencer = rigidbody::sequencer::SequenceParser().parse_text(script_obj->script);
    sequencer->run();
    return 0;
}, status);}