#pragma once

#include <md/programs/options/forcefields/IForcefield.h>

namespace ausaxs::md::option::force_field {
    struct AMBER99SB : IForcefield {
        std::string filename() const override { return "amber99sb"; }
    };
}