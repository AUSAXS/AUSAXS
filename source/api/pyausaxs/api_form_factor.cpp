// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/pyausaxs/api_form_factor.h>
#include <api/ObjectStorage.h>
#include <form_factor/FormFactorType.h>
#include <form_factor/FormFactorTable.h>
#include <form_factor/ExvTable.h>

#include <functional>
#include <array>

using namespace ausaxs;
using namespace ausaxs::form_factor;

namespace {
    std::tuple<std::reference_wrapper<const std::array<double, 5>>, std::reference_wrapper<const std::array<double, 5>>, std::reference_wrapper<const double>> 
        get_gaussian_coefficients(form_factor_t type)
    {
        switch (type) {
            case form_factor_t::H:      return {std::cref(constants::form_factor::H::a),        std::cref(constants::form_factor::H::b),        std::cref(constants::form_factor::H::c)};
            case form_factor_t::C:      return {std::cref(constants::form_factor::C::a),        std::cref(constants::form_factor::C::b),        std::cref(constants::form_factor::C::c)};
            case form_factor_t::N:      return {std::cref(constants::form_factor::N::a),        std::cref(constants::form_factor::N::b),        std::cref(constants::form_factor::N::c)};
            case form_factor_t::O:      return {std::cref(constants::form_factor::O::a),        std::cref(constants::form_factor::O::b),        std::cref(constants::form_factor::O::c)};
            case form_factor_t::S:      return {std::cref(constants::form_factor::S::a),        std::cref(constants::form_factor::S::b),        std::cref(constants::form_factor::S::c)};
            case form_factor_t::CH:     return {std::cref(constants::form_factor::CH_sp3::a),   std::cref(constants::form_factor::CH_sp3::b),   std::cref(constants::form_factor::CH_sp3::c)};
            case form_factor_t::CH2:    return {std::cref(constants::form_factor::CH2_sp3::a),  std::cref(constants::form_factor::CH2_sp3::b),  std::cref(constants::form_factor::CH2_sp3::c)};
            case form_factor_t::CH3:    return {std::cref(constants::form_factor::CH3_sp3::a),  std::cref(constants::form_factor::CH3_sp3::b),  std::cref(constants::form_factor::CH3_sp3::c)};
            case form_factor_t::NH:     return {std::cref(constants::form_factor::NH::a),       std::cref(constants::form_factor::NH::b),       std::cref(constants::form_factor::NH::c)};
            case form_factor_t::NH2:    return {std::cref(constants::form_factor::NH2::a),      std::cref(constants::form_factor::NH2::b),      std::cref(constants::form_factor::NH2::c)};
            case form_factor_t::NH3:    return {std::cref(constants::form_factor::NH3_plus::a), std::cref(constants::form_factor::NH3_plus::b), std::cref(constants::form_factor::NH3_plus::c)};
            case form_factor_t::OH:     return {std::cref(constants::form_factor::OH_alc::a),   std::cref(constants::form_factor::OH_alc::b),   std::cref(constants::form_factor::OH_alc::c)};
            case form_factor_t::SH:     return {std::cref(constants::form_factor::SH::a),       std::cref(constants::form_factor::SH::b),       std::cref(constants::form_factor::SH::c)};
            case form_factor_t::OTHER:  return {std::cref(constants::form_factor::other::a),    std::cref(constants::form_factor::other::b),    std::cref(constants::form_factor::other::c)};
            case form_factor_t::EXCLUDED_VOLUME: return {std::cref(constants::form_factor::excluded_volume::a), std::cref(constants::form_factor::excluded_volume::b), std::cref(constants::form_factor::excluded_volume::c)};
            default: throw std::runtime_error("get_gaussian_coefficients: Unknown form factor type (enum " + std::to_string(static_cast<int>(type)) + ")");
        }
    } 
}

struct _ff_valid_form_factor_types_obj {
    std::vector<std::string> types;
    std::vector<const char*> types_ptr;
};
int ff_valid_form_factor_types(
    const char*** types,
    int* n_types,
    int* status
) {return execute_with_catch([&]() {
    _ff_valid_form_factor_types_obj obj;
    for (unsigned int i = 0; i < form_factor::get_count(); ++i) {
        obj.types.emplace_back(form_factor::to_string(static_cast<form_factor::form_factor_t>(i)));
    }
    obj.types_ptr.resize(obj.types.size());
    for (size_t i = 0; i < obj.types.size(); ++i) {
        obj.types_ptr[i] = obj.types[i].c_str();
    }
    *n_types = static_cast<int>(obj.types.size());
    int id = api::ObjectStorage::register_object(std::move(obj));
    auto ref = api::ObjectStorage::get_object<_ff_valid_form_factor_types_obj>(id);
    *types = ref->types_ptr.data();
    *status = 0;
    return id;
}, status);}

void ff_get_five_gaussian_coefficients(
    const char* element, 
    double* a, double* b, double* c,
    int* status
) {return execute_with_catch([&]() {
    form_factor::form_factor_t type = from_string(element);
    auto[a_ref, b_ref, c_ref] = get_gaussian_coefficients(type);
    std::copy(a_ref.get().begin(), a_ref.get().end(), a);
    std::copy(b_ref.get().begin(), b_ref.get().end(), b);
    *c = c_ref.get();
}, status);}

void ff_get_current_exv_volume(
    const char* element,
    double* volume,
    int* status
) {return execute_with_catch([&]() {
    form_factor::form_factor_t type = from_string(element);
    const auto& exv_table = constants::exv::get_exv_set();
    switch (type) {
        case form_factor::form_factor_t::H: *volume = exv_table.H; break;
        case form_factor::form_factor_t::C: *volume = exv_table.C; break;
        case form_factor::form_factor_t::N: *volume = exv_table.N; break;
        case form_factor::form_factor_t::O: *volume = exv_table.O; break;
        case form_factor::form_factor_t::S: *volume = exv_table.S; break;
        case form_factor::form_factor_t::CH: *volume = exv_table.CH; break;
        case form_factor::form_factor_t::CH2: *volume = exv_table.CH2; break;
        case form_factor::form_factor_t::CH3: *volume = exv_table.CH3; break;
        case form_factor::form_factor_t::NH: *volume = exv_table.NH; break;
        case form_factor::form_factor_t::NH2: *volume = exv_table.NH2; break;
        case form_factor::form_factor_t::NH3: *volume = exv_table.NH3; break;
        case form_factor::form_factor_t::OH: *volume = exv_table.OH; break;
        case form_factor::form_factor_t::SH: *volume = exv_table.SH; break;
        case form_factor::form_factor_t::OTHER: *volume = constants::exv::Ar; break;
        default: throw std::runtime_error("ff_get_current_exv_volume: Unknown form factor type \"" + std::string(element) + "\"");
    }
}, status);}
