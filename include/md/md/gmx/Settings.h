#pragma once

#include <gmx/SmartOption.h>

namespace gmx {
    struct setting {
        inline static SmartOption<std::string> gmx_path = {"gmx", {"gmx_exe", "gmx_executable", "gmx"}};
        inline static SmartOption<std::string> buffer_path = {"", {"buffer_path", "buffer"}};
        inline static SmartOption<std::string> water_model = {"tip4p2005", {"water_model", "water"}};
        inline static SmartOption<std::string> force_field = {"amber99sb-ildn", {"force_field", "forcefield", "ff"}};
        inline static SmartOption<std::string> minimization_sim_location = {"lucy", {"minimization_sim_location", "minimization_sim", "minimization"}};
        inline static SmartOption<std::string> box_type = {"dodecahedron", {"box_type", "boxtype", "box"}};
        inline static SmartOption<std::string> cation = {"Na", {"cation", "cation_type", "cationtype"}};
        inline static SmartOption<std::string> anion = {"Cl", {"anion", "anion_type", "aniontype"}};

        inline static SmartOption<std::string> thermalization_sim_location = {"smaug", {"thermalization_sim_location", "thermalization_sim", "thermalization"}};
        inline static SmartOption<std::string> production_sim_location = {"smaug", {"production_sim_location", "production_sim", "production"}};
    };
}