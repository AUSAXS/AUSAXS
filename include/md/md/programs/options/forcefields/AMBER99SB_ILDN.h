#pragma once

#include <md/programs/options/forcefields/IForcefield.h>

namespace ausaxs::md::option::force_field {
    struct AMBER99SB_ILDN : IForcefield {
        std::string filename() const override { return "amber99sb-ildn"; }
    };
}