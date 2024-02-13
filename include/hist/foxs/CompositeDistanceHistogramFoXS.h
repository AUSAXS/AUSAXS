#pragma once

#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvg.h>
#include <hist/distribution/Distribution3D.h>
#include <hist/distribution/Distribution2D.h>
#include <hist/distribution/Distribution1D.h>
#include <hist/distribution/WeightedDistribution3D.h>
#include <hist/distribution/WeightedDistribution2D.h>
#include <hist/distribution/WeightedDistribution1D.h>

namespace hist {
    /**
     * @brief An alternative to CompositeDistanceHistogramFFExplicit that uses the same form factors as FoXS. 
     */
    class CompositeDistanceHistogramFoXS : public CompositeDistanceHistogramFFAvg {
        public: 
            /**
             * @brief Create an unweighted composite distance histogram with form factors.
             * 
             * @param p_aa The partial distance histogram for atom-atom interactions
             * @param p_ax The partial distance histogram for atom-excluded volume interactions
             * @param p_xx The partial distance histogram for excluded volume-excluded volume interactions
             * @param p_aw The partial distance histogram for atom-water interactions
             * @param p_wx The partial distance histogram for water-excluded volume interactions
             * @param p_ww The partial distance histogram for water-water interactions
             * @param p_tot The total distance histogram. This is only used for determining the maximum distance.
             */
            CompositeDistanceHistogramFoXS(
                hist::Distribution3D&& p_aa, 
                hist::Distribution3D&& p_ax, 
                hist::Distribution3D&& p_xx, 
                hist::Distribution2D&& p_aw, 
                hist::Distribution2D&& p_wx, 
                hist::Distribution1D&& p_ww,
                hist::Distribution1D&& p_tot
            );

            /**
             * @brief Create a weighted composite distance histogram with form factors.
             * 
             * @param p_aa The partial distance histogram for atom-atom interactions
             * @param p_ax The partial distance histogram for atom-excluded volume interactions
             * @param p_xx The partial distance histogram for excluded volume-excluded volume interactions
             * @param p_aw The partial distance histogram for atom-water interactions
             * @param p_wx The partial distance histogram for water-excluded volume interactions
             * @param p_ww The partial distance histogram for water-water interactions
             * @param p_tot The total distance histogram. This is only used to extract the bin centers.
             */
            CompositeDistanceHistogramFoXS(
                hist::Distribution3D&& p_aa, 
                hist::Distribution3D&& p_ax, 
                hist::Distribution3D&& p_xx, 
                hist::Distribution2D&& p_aw, 
                hist::Distribution2D&& p_wx, 
                hist::Distribution1D&& p_ww, 
                hist::WeightedDistribution1D&& p_tot
            );

            ~CompositeDistanceHistogramFoXS() override;

            ScatteringProfile debye_transform() const override; // @copydoc DistanceHistogram::debye_transform() const

            /**
             * @brief Get the intensity profile for atom-atom interactions.
             */
            virtual ScatteringProfile get_profile_ax() const override;

            /**
             * @brief Get the intensity profile for atom-water interactions.
             */
            virtual ScatteringProfile get_profile_xx() const override;

            /**
             * @brief Get the intensity profile for water-water interactions.
             */
            virtual ScatteringProfile get_profile_wx() const override;

            /**
             * @brief Get the intensity profile for atom-water interactions.
             */
            virtual ScatteringProfile get_profile_aw() const override;

            /**
             * @brief Get the intensity profile for water-water interactions.
             */
            virtual ScatteringProfile get_profile_ww() const override;

            /**
             * @brief Get the intensity profile for atom-atom interactions.
             */
            virtual ScatteringProfile get_profile_aa() const override;

        protected:
            double G_factor(double q) const;

        private:
            container::Container3D<double> cp_ax;
            container::Container3D<double> cp_xx;
            container::Container2D<double> cp_wx;
    };
}