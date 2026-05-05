// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/pyausaxs/api_settings.h>
#include <api/ObjectStorage.h>
#include <settings/SettingsIO.h>
#include <settings/SettingRef.h>
#include <settings/All.h>

using namespace ausaxs;

struct _get_setting_obj {
    std::string value;
    std::string type;
};
int get_setting(
    const char* name,
    const char** value,
    const char** type,
    int* status
) {return execute_with_catch([&]() {
    std::string name_str(name);
    const auto& map = settings::io::detail::ISettingRef::get_stored_settings();
    if (!map.contains(name_str)) {ErrorMessage::last_error = "Unknown setting: \"" + name_str + "\""; return -1;}
    const auto& setting = map.at(name_str);
    auto obj_id = api::ObjectStorage::register_object(_get_setting_obj{
        .value = setting->get(),
        .type = setting->type()
    });
    auto obj = api::ObjectStorage::get_object<_get_setting_obj>(obj_id);
    *value = obj->value.c_str();
    *type = obj->type.c_str();
    return obj_id;
}, status);}

void set_setting(
    const char* name,
    const char* value,
    int* status
) {return execute_with_catch([&]() {
    const auto& map = settings::io::detail::ISettingRef::get_stored_settings();
    if (!map.contains(name)) {ErrorMessage::last_error = "Unknown setting: \"" + std::string(name) + "\""; return;}
    map.at(name)->set({value});
}, status);}