#pragma once

#include <utility/Settings.h>
#include <utility/Exceptions.h>

namespace setting {
    void set_miller_generation_strategy(const std::string& strategy);

    void set_rigidbody_parameter_strategy(const std::string& strategy);

    void set_rigidbody_transformation_strategy(const std::string& strategy);

    void set_rigidbody_selection_strategy(const std::string& strategy);

    void set_hydration_placement_strategy(const std::string& strategy) {
        if (strategy == "Radial") {setting::grid::placement_strategy = setting::grid::PlacementStrategy::RadialStrategy;}
        else if (strategy == "Axes") {setting::grid::placement_strategy = setting::grid::PlacementStrategy::AxesStrategy;}
        else if (strategy == "Jan") {setting::grid::placement_strategy = setting::grid::PlacementStrategy::JanStrategy;}
        else {throw except::invalid_argument("StrategySettings::Unknown strategy: \"" + strategy + "\"");}
    }

    void set_hydration_culling_strategy(const std::string& strategy);
}