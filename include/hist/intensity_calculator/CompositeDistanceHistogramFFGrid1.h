#include <hist/intensity_calculator/CompositeDistanceHistogramFFGrid.h>

namespace hist {
    class CompositeDistanceHistogramFFGrid1 : public CompositeDistanceHistogramFFGrid {
        public:
            using CompositeDistanceHistogramFFGrid::CompositeDistanceHistogramFFGrid;
            ~CompositeDistanceHistogramFFGrid1() override = default;

        private:
            virtual double exv_factor(double) const override {
                return cx;
            }
    };
}