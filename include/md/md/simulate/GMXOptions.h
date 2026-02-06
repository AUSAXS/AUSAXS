#pragma once

#include <md/programs/gmx.h>
#include <md/programs/options/forcefields/IForcefield.h>
#include <md/programs/options/water_models/IWaterModel.h>

#include <unordered_map>

namespace ausaxs::md {
    struct SystemSettings {
        option::Forcefield forcefield;
        option::WaterModel watermodel;
        option::BoxType boxtype;
        option::Cation cation;
        option::Anion anion;
        std::unordered_map<std::string, std::string> custom_options;
    };
}