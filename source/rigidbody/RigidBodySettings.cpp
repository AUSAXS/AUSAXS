#include <rigidbody/RigidBodySettings.h>

namespace settings::rigidbody {
    settings::detail::SmartOption<unsigned int> iterations(1000, "rigidbody-iterations");
    settings::detail::SmartOption<double> bond_distance(3, "rigidbody-bond-distance");

    namespace detail {
        settings::detail::SmartOption<std::vector<int>> constraints({}, "rigidbody-constraints");
        settings::detail::SmartOption<std::string> calibration_file("", {"calibration", "rigidbody-calibration-file"});
    }
}