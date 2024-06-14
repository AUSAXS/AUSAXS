#pragma once

#include <hist/distribution/DistributionFwd.h>

#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvg.h>

namespace hist {
    /**
     * @brief A class containing partial distance histograms for the different types of interactions and atomic types. 
     *        Beyond the functionality of CompositeDistanceHistogram, this class also uses individual form factors for each atomic type.
     *        Unique exluded volume form factors are used for each atomic type.
     *        For more information, see CompositeDistanceHistogram.
     */
    template<typename AAFormFactorTableType, typename AXFormFactorTableType, typename XXFormFactorTableType>
    class CompositeDistanceHistogramFFExplicitBase : public CompositeDistanceHistogramFFAvgBase<AAFormFactorTableType> {
        public: 
            CompositeDistanceHistogramFFExplicitBase();
            CompositeDistanceHistogramFFExplicitBase(CompositeDistanceHistogramFFExplicitBase&&) noexcept;
            CompositeDistanceHistogramFFExplicitBase& operator=(CompositeDistanceHistogramFFExplicitBase&&) noexcept;
            virtual ~CompositeDistanceHistogramFFExplicitBase() override;

            /**
             * @brief Create a new unweighted composite distance histogram with form factors.
             * 
             * @param p_aa The partial distance histogram for atom-atom interactions.
             * @param p_ax The partial distance histogram for atom-excluded volume interactions.
             * @param p_xx The partial distance histogram for excluded volume-excluded volume interactions.
             * @param p_aw The partial distance histogram for atom-water interactions.
             * @param p_wx The partial distance histogram for water-excluded volume interactions.
             * @param p_ww The partial distance histogram for water-water interactions.
             * @param p_tot The total distance histogram. This is only used for determining the maximum distance.
             */
            CompositeDistanceHistogramFFExplicitBase(
                hist::Distribution3D&& p_aa, 
                hist::Distribution3D&& p_ax, 
                hist::Distribution3D&& p_xx, 
                hist::Distribution2D&& p_aw, 
                hist::Distribution2D&& p_wx, 
                hist::Distribution1D&& p_ww,
                hist::Distribution1D&& p_tot
            );

            /**
             * @brief Create a new weighted composite distance histogram with form factors.
             * 
             * @param p_aa The partial distance histogram for atom-atom interactions.
             * @param p_ax The partial distance histogram for atom-excluded volume interactions.
             * @param p_xx The partial distance histogram for excluded volume-excluded volume interactions.
             * @param p_aw The partial distance histogram for atom-water interactions.
             * @param p_wx The partial distance histogram for water-excluded volume interactions.
             * @param p_ww The partial distance histogram for water-water interactions.
             * @param p_tot The total distance histogram. This is only used to extract the bin centers.
             */
            CompositeDistanceHistogramFFExplicitBase(
                hist::Distribution3D&& p_aa, 
                hist::Distribution3D&& p_ax, 
                hist::Distribution3D&& p_xx, 
                hist::Distribution2D&& p_aw, 
                hist::Distribution2D&& p_wx, 
                hist::Distribution1D&& p_ww, 
                hist::WeightedDistribution1D&& p_tot
            );

            // @copydoc DistanceHistogram::debye_transform() const
            ScatteringProfile debye_transform() const override;

            // @copydoc DistanceHistogram::debye_transform(const std::vector<double>&) const
            virtual SimpleDataset debye_transform(const std::vector<double>& q) const override;

            virtual ScatteringProfile get_profile_ax() const override; // @copydoc ICompositeDistanceHistogramExv::get_profile_ax() const
            virtual ScatteringProfile get_profile_xx() const override; // @copydoc ICompositeDistanceHistogramExv::get_profile_xx() const
            virtual ScatteringProfile get_profile_wx() const override; // @copydoc ICompositeDistanceHistogramExv::get_profile_wx() const

            const AAFormFactorTableType get_ffaa_table() const;
            virtual const AXFormFactorTableType& get_ffax_table() const = 0;
            virtual const XXFormFactorTableType& get_ffxx_table() const = 0;

        protected:
            // @copydoc CompositeDistanceHistogramFFAvgBase::exv_factor(double) const
            double exv_factor(double q) const override;

            struct {hist::Distribution3D xx, ax; hist::Distribution2D wx;} exv_distance_profiles;
    };
}