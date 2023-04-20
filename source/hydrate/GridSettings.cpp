#include <hydrate/GridSettings.h>

using namespace settings::detail;

namespace settings::grid {
    SmartOption<double> percent_water(0.1, "percent_water");
    SmartOption<double> ra(2.4, {"atomic-radius", "ra"});
    SmartOption<double> rh(1.5, {"water-radius", "rh"});
    SmartOption<double> ra_effective(2.4, "effective-radius");
    SmartOption<double> width(1, {"grid-width"});
    SmartOption<double> scaling(0.25, "grid-scaling");
    SmartOption<bool> cubic(false, "grid-cubic");
    SmartOption<Limit3D> axes(Limit3D(-250, 250, -250, 250, -250, 250), "grid-axes");

    namespace detail {
        SmartOption<double> min_score(0.1, "min-score");
    }
}