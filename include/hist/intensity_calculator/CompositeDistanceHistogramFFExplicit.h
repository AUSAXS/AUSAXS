#pragma once

#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvg.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/distribution/GenericDistribution2D.h>
#include <hist/distribution/GenericDistribution3D.h>

#include <vector>

namespace hist {
    /**
     * @brief A class containing multiple partial distance histograms for multiple form factors. 
     */
    class CompositeDistanceHistogramFFExplicit : public CompositeDistanceHistogramFFAvg {
        public: 
            CompositeDistanceHistogramFFExplicit();

            /**
             * @brief Construct a new Composite Distance Histogram FF object
             * 
             * @param p_aa Partial distance histogram for atom-atom interactions
             * @param p_ax Partial distance histogram for atom-excluded volume interactions
             * @param p_xx Partial distance histogram for excluded volume-excluded volume interactions
             * @param p_aw Partial distance histogram for atom-water interactions
             * @param p_wx Partial distance histogram for water-excluded volume interactions
             * @param p_ww Partial distance histogram for water-water interactions
             * @param p_tot Total distance histogram
             * @param axis Distance axis
             */
            CompositeDistanceHistogramFFExplicit(
                hist::Distribution3D&& p_aa, 
                hist::Distribution3D&& p_ax, 
                hist::Distribution3D&& p_xx, 
                hist::Distribution2D&& p_aw, 
                hist::Distribution2D&& p_wx, 
                hist::Distribution1D&& p_ww, 
                const Axis& axis
            );

            CompositeDistanceHistogramFFExplicit(
                hist::WeightedDistribution3D&& p_aa, 
                hist::WeightedDistribution3D&& p_ax, 
                hist::WeightedDistribution3D&& p_xx, 
                hist::WeightedDistribution2D&& p_aw, 
                hist::WeightedDistribution2D&& p_wx, 
                hist::WeightedDistribution1D&& p_ww, 
                const Axis& axis
            );

            ~CompositeDistanceHistogramFFExplicit() override;

            ScatteringProfile debye_transform() const override;

            /**
             * @brief Get the intensity profile for atom-atom interactions.
             */
            virtual const ScatteringProfile get_profile_ax() const;

            /**
             * @brief Get the intensity profile for atom-water interactions.
             */
            virtual const ScatteringProfile get_profile_xx() const;

            /**
             * @brief Get the intensity profile for water-water interactions.
             */
            virtual const ScatteringProfile get_profile_wx() const;

        protected:
            double exv_factor(double q) const override;

        private:
            hist::Distribution3D cp_ax;
            hist::Distribution3D cp_xx;
            hist::Distribution2D cp_wx;
    };
}