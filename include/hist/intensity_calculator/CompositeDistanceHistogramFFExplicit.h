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
    template<bool use_weighted_distribution>
    class CompositeDistanceHistogramFFExplicit : public CompositeDistanceHistogramFFAvg<use_weighted_distribution> {
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
                hist::GenericDistribution3D<use_weighted_distribution>&& p_aa, 
                hist::GenericDistribution3D<use_weighted_distribution>&& p_ax, 
                hist::GenericDistribution3D<use_weighted_distribution>&& p_xx, 
                hist::GenericDistribution2D<use_weighted_distribution>&& p_wa, 
                hist::GenericDistribution2D<use_weighted_distribution>&& p_wx, 
                hist::GenericDistribution1D<use_weighted_distribution>&& p_ww, 
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
            double G_factor(double q) const;

        private:
            hist::GenericDistribution3D<use_weighted_distribution> cp_ax;
            hist::GenericDistribution3D<use_weighted_distribution> cp_xx;
            hist::GenericDistribution2D<use_weighted_distribution> cp_wx;
    };
}