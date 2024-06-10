#include <md/programs/gmx.h>

using namespace gmx;

std::string option::to_string(Forcefield opt) {
    switch (opt) {
        case Forcefield::AMBER99SB:
            return "amber99sb";
        case Forcefield::AMBER99SB_ILDN:
            return "amber99sb-ildn";
        default:
            throw except::unknown_type("gmx::to_string: Unknown forcefield. (Did you forget to add it to the enum?)");
    }
}

std::string option::to_string(WaterModel opt) {
    switch (opt) {
        case WaterModel::TIP3P:
            return "tip3p";
        case WaterModel::TIP4P:
            return "tip4p";
        case WaterModel::TIP4P2005:
            return "tip4p2005";
        default:
            throw except::unknown_type("gmx::to_string: Unknown water model. (Did you forget to add it to the enum?)");
    }
}

std::string option::to_string(BoxType opt) {
    switch (opt) {
        case BoxType::CUBIC:
            return "cubic";
        case BoxType::TRICLINIC:
            return "triclinic";
        case BoxType::DODECAHEDRON:
            return "dodecahedron";
        case BoxType::OCTAHEDRON:
            return "octahedron";
        default:
            throw except::unknown_type("gmx::to_string: Unknown box type. (Did you forget to add it to the enum?)");
    }
}

std::string option::to_string(Cation opt) {
    switch (opt) {
        case Cation::NA:
            return "NA";
        default:
            throw except::unknown_type("gmx::to_string: Unknown cation type. (Did you forget to add it to the enum?)");
    }
}

std::string option::to_string(Anion opt) {
    switch (opt) {
        case Anion::CL:
            return "CL";
        default:
            throw except::unknown_type("gmx::to_string: Unknown anion type. (Did you forget to add it to the enum?)");
    }
}