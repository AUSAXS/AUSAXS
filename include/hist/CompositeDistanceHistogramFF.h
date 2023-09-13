#pragma once

#include <hist/DistanceHistogram.h>

#include <vector>

namespace hist {
    class Histogram;
    class CompositeDistanceHistogram;

    /**
     * @brief A class containing multiple partial distance histograms for multiple form factors. 
     */
    class CompositeDistanceHistogramFF : public DistanceHistogram {
        public: 
            CompositeDistanceHistogramFF() = default;

            ~CompositeDistanceHistogramFF() override;

            Histogram debye_transform() const override;

            Histogram debye_transform(const std::vector<double>& q) const override;

            void apply_water_scaling_factor(double k);

        private:
            std::vector<CompositeDistanceHistogram> cdhs;
    };
}