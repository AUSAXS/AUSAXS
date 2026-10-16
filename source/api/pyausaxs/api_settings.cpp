// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/pyausaxs/api_settings.h>

#include <api/ObjectStorage.h>
#include <constants/ConstantsAxes.h>
#include <settings/All.h>
#include <settings/SettingRef.h>
#include <settings/SettingsIO.h>

using namespace ausaxs;

namespace {
    struct _get_setting_obj {
        std::string value;
        std::string type;
    };
}
int get_setting(
    const char* name,
    const char** value,
    const char** type,
    int* status
) {return execute_with_catch([&]() {
    std::string name_str(name);
    const auto& map = settings::io::detail::ISettingRef::get_stored_settings();
    if (!map.contains(name_str)) {throw except::invalid_argument("Unknown setting: \"" + name_str + "\"");}
    const auto& setting = map.at(name_str);
    auto obj_id = api::ObjectStorage::register_object(_get_setting_obj{
        .value = setting->get(),
        .type = setting->type()
    });
    auto* obj = api::ObjectStorage::get_object<_get_setting_obj>(obj_id);
    *value = obj->value.c_str();
    *type = obj->type.c_str();
    return obj_id;
}, status);}

void set_setting(
    const char* name,
    const char* value,
    int* status
) {execute_with_catch([&]() {
    const auto& map = settings::io::detail::ISettingRef::get_stored_settings();
    if (!map.contains(name)) {throw except::invalid_argument("Unknown setting: \"" + std::string(name) + "\"");}
    map.at(name)->set({value});
}, status);}

void save_settings(
    const char* path,
    int* status
) {execute_with_catch([&]() {
    settings::write(io::File(path));
}, status);}

void load_settings(
    const char* path,
    int* status
) {execute_with_catch([&]() {
    settings::read(io::ExistingFile(path));
}, status);}

namespace {
    struct _get_q_axis_obj {
        std::vector<double> q;
    };
}
int get_q_axis(
    const double** q,
    int* n_points,
    int* status
) {return execute_with_catch([&]() {
    auto obj_id = api::ObjectStorage::register_object(_get_q_axis_obj{
        .q = constants::axes::q_axis.sub_axis_covering(settings::axes::qmin, settings::axes::qmax).as_vector()
    });
    auto* obj = api::ObjectStorage::get_object<_get_q_axis_obj>(obj_id);
    *q = obj->q.data();
    *n_points = static_cast<int>(obj->q.size());
    return obj_id;
}, status);}
