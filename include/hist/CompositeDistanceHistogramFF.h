#pragma once

#include <hist/CompositeDistanceHistogram.h>
#include <utility/Container1D.h>
#include <utility/Container2D.h>
#include <utility/Container3D.h>

#include <vector>

namespace hist {
    class Histogram;

    /**
     * @brief A class containing multiple partial distance histograms for multiple form factors. 
     */
    class CompositeDistanceHistogramFF : public CompositeDistanceHistogram {
        public: 
            CompositeDistanceHistogramFF() = default;

            CompositeDistanceHistogramFF(Container3D<double>&& p_pp, Container2D<double>&& p_hp, Container1D<double>&& p_hh, std::vector<double>&& p_tot, const Axis& axis);

            ~CompositeDistanceHistogramFF() override;

            Histogram debye_transform() const override;

            Histogram debye_transform(const std::vector<double>& q) const override;

            void apply_water_scaling_factor(double k);

        private:
            Container3D<double> p_pp;
            Container2D<double> p_hp;
            Container1D<double> p_hh;
    };
}