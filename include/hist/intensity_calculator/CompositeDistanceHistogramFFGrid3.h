#include <hist/intensity_calculator/CompositeDistanceHistogramFFGrid.h>

namespace hist {
    class CompositeDistanceHistogramFFGrid3 : public CompositeDistanceHistogramFFGrid {
        public:
            using CompositeDistanceHistogramFFGrid::CompositeDistanceHistogramFFGrid;
            ~CompositeDistanceHistogramFFGrid3() override = default;

        private:
            virtual double exv_factor(double q) const override {
                // G(q) factor from CRYSOL
                constexpr double rm = 1.62;
                constexpr double c = constexpr_math::pow(4*constants::pi/3, 3./2)*constants::pi*rm*rm*constants::form_factor::s_to_q_factor;
                return std::clamp(std::pow(cx, 3)*std::exp(-c*(cx*cx - 1)*q*q), 0., 10.);
            }
    };
}