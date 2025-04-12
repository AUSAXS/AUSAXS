#pragma once

#include <md/programs/gmx.h>
#include <md/utility/files/all.h>
#include <md/programs/mdrun/Execution.h>
#include <md/utility/files/MDPCreator.h>
#include <md/programs/options/forcefields/IForcefield.h>
#include <md/programs/options/water_models/IWaterModel.h>

namespace ausaxs::md {
    struct SystemSettings {
        option::Forcefield forcefield;
        option::WaterModel watermodel;
        option::BoxType boxtype;
        option::Cation cation;
        option::Anion anion;
    };
}