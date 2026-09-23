// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/pyausaxs/api_form_factor.h>

#include <api/ObjectStorage.h>
#include <form_factor/FormFactorType.h>
#include <form_factor/lookup/ExvTableManager.h>

#include <algorithm>

using namespace ausaxs;
using namespace ausaxs::form_factor;

namespace {
struct _ff_valid_form_factor_types_obj {
    std::vector<std::string> types;
    std::vector<const char*> types_ptr;
};
}
int ff_valid_form_factor_types(
    const char*** types,
    int* n_types,
    int* status
) {return execute_with_catch([&]() {
    _ff_valid_form_factor_types_obj obj;
    for (int i = 0; i < form_factor::total_ff_count; ++i) {
        obj.types.emplace_back(form_factor::to_string(static_cast<form_factor::form_factor_t>(i)));
    }
    obj.types_ptr.resize(obj.types.size());
    for (size_t i = 0; i < obj.types.size(); ++i) {
        obj.types_ptr[i] = obj.types[i].c_str();
    }
    *n_types = static_cast<int>(obj.types.size());
    int id = api::ObjectStorage::register_object(std::move(obj));
    auto* ref = api::ObjectStorage::get_object<_ff_valid_form_factor_types_obj>(id);
    *types = ref->types_ptr.data();
    *status = 0;
    return id;
}, status);}

void ff_get_five_gaussian_coefficients(
    const char* element, 
    double* a, double* b, double* c,
    int* status
) {execute_with_catch([&]() {
    form_factor::form_factor_t type = from_string(element);
    const auto& coefficients = form_factor::get_info(type).coefficients;
    std::ranges::copy(coefficients.a, a);
    std::ranges::copy(coefficients.b, b);
    *c = coefficients.c;
}, status);}

void ff_get_current_exv_volume(
    const char* element,
    double* volume,
    int* status
) {execute_with_catch([&]() {
    form_factor::form_factor_t type = from_string(element);
    *volume = ExvTableManager::get_current_exv_table()->get(type);
}, status);}
