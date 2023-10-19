#pragma once

#include <hist/CompositeDistanceHistogramFF.h>
#include <container/Container1D.h>
#include <container/Container2D.h>
#include <container/Container3D.h>

#include <vector>

namespace hist {
    /**
     * @brief A class containing multiple partial distance histograms for multiple form factors. 
     */
    class CompositeDistanceHistogramFFExplicit : public CompositeDistanceHistogramFF {
        public: 
            CompositeDistanceHistogramFFExplicit();

            /**
             * @brief Construct a new Composite Distance Histogram FF object
             * 
             * @param p_pp Partial distance histogram for atom-atom interactions
             * @param p_hp Partial distance histogram for atom-water interactions
             * @param p_hh Partial distance histogram for water-water interactions
             * @param p_tot Total distance histogram
             * @param axis Distance axis
             * @param Z_exv_avg Average charge of excluded volume
             */
            CompositeDistanceHistogramFFExplicit(container::Container3D<double>&& p_aa, container::Container3D<double>&& p_ax, container::Container3D<double>&& p_xx, 
                                                 container::Container2D<double>&& p_wa, container::Container2D<double>&& p_wx, container::Container1D<double>&& p_ww, 
                                                 std::vector<double>&& p_tot, const Axis& axis);

            ~CompositeDistanceHistogramFFExplicit() override;

            ScatteringProfile debye_transform() const override;

        private:
            container::Container3D<double> cp_ax;
            container::Container3D<double> cp_xx;
            container::Container2D<double> cp_wx;
    };
}