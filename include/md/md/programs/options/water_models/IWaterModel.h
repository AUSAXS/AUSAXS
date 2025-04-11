#pragma once

#include <md/programs/options/forcefields/IForcefield.h>
#include <utility/observer_ptr.h>
#include <io/IOFwd.h>

#include <string>
#include <memory>

namespace ausaxs::md::option {
    enum class WaterModel {
        TIP3P,
        TIP4P,
        TIP4PEw,
        TIP4P2005,
        TIP5P,
        SPC,
        SPCE,
    };

    struct IWaterModel {
        static std::unique_ptr<IWaterModel> construct(WaterModel wm);

        virtual std::string filename() const = 0;
        virtual std::string name() const = 0;
        virtual std::string info() const = 0;
        virtual std::string get_file_content() const = 0;
        bool exists(observer_ptr<option::IForcefield> ff) const;
        void create(observer_ptr<option::IForcefield> ff) const;
    };
}