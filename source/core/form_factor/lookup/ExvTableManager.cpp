#pragma once

#include <form_factor/lookup/ExvTableManager.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <settings/ExvSettings.h>

#include <unordered_map>

using namespace ausaxs;
using namespace ausaxs::form_factor; 

observer_ptr<const constants::exv::detail::ExvSet> ExvTableManager::get_current_exv_table() {
    switch (settings::exv::exv_set) {
        case settings::exv::ExvSet::Traube: return &constants::exv::Traube;
        case settings::exv::ExvSet::Voronoi_explicit_H: return &constants::exv::Voronoi_explicit_H;
        case settings::exv::ExvSet::Voronoi_implicit_H: return &constants::exv::Voronoi_implicit_H;
        case settings::exv::ExvSet::MinimumFluctutation_explicit_H: return &constants::exv::MinimumFluctuation_explicit_H;
        case settings::exv::ExvSet::MinimumFluctutation_implicit_H: return &constants::exv::MinimumFluctuation_implicit_H;
        case settings::exv::ExvSet::vdw: return &constants::exv::vdw;
        case settings::exv::ExvSet::Custom: {
            assert(custom_exv_tables && "Custom excluded volume tables must be set before they can be used.");
            return custom_exv_tables.get();
        }
        default: 
            throw std::runtime_error(
                "constants::displaced_volume::get_exv_set: Invalid displaced volume set" 
                "(enum " + std::to_string(static_cast<int>(settings::exv::exv_set)) + ")");
    }
}

const detail::ExvFormFactorSet& ExvTableManager::get_current_exv_form_factor_set() {
    static auto available_sets = std::unordered_map<settings::exv::ExvSet, detail::ExvFormFactorSet>{
        {settings::exv::ExvSet::Default, lookup::exv::standard}
    };

    // always update custom sets to ensure they reflect changes to the exv table
    if (settings::exv::exv_set == settings::exv::ExvSet::Custom) {
        //? add caching for custom tables? 
        available_sets[settings::exv::ExvSet::Custom] = detail::ExvFormFactorSet(*get_current_exv_table());
        return available_sets[settings::exv::ExvSet::Custom];
    } else {
        // check if the current set is already available, otherwise create it
        if (auto it = available_sets.find(settings::exv::exv_set); it != available_sets.end()) {
            return it->second;
        }
        available_sets[settings::exv::exv_set] = detail::ExvFormFactorSet(*get_current_exv_table());
        return available_sets[settings::exv::exv_set];
    }
}

void ExvTableManager::set_custom_exv_table(const constants::exv::detail::ExvSet& set) {
    custom_exv_tables = std::make_unique<constants::exv::detail::ExvSet>(set);
    FormFactorManager::_needs_refresh = true;
}