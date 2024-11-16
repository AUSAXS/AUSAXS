#pragma once

#include <hist/intensity_calculator/CompositeDistanceHistogramFFExplicit.h>
#include <hist/intensity_calculator/crysol/FormFactorCrysol.h>

namespace ausaxs::hist {
    /**
     * @brief An alternative to CompositeDistanceHistogramFFExplicit that mimics the CRYSOL excluded volume fitting. 
     */
    class CompositeDistanceHistogramCrysol : public CompositeDistanceHistogramFFExplicitBase<form_factor::storage::atomic::table_t, form_factor::storage::cross::table_t, form_factor::storage::exv::table_t>{
        public:
            CompositeDistanceHistogramCrysol();
            CompositeDistanceHistogramCrysol(const CompositeDistanceHistogramCrysol&);
            CompositeDistanceHistogramCrysol(CompositeDistanceHistogramCrysol&&) noexcept;
            CompositeDistanceHistogramCrysol& operator=(CompositeDistanceHistogramCrysol&&) noexcept;
            CompositeDistanceHistogramCrysol& operator=(const CompositeDistanceHistogramCrysol&);
            virtual ~CompositeDistanceHistogramCrysol() override;

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
            CompositeDistanceHistogramCrysol(
                hist::Distribution3D&& p_aa, 
                hist::Distribution3D&& p_ax, 
                hist::Distribution3D&& p_xx, 
                hist::Distribution2D&& p_aw, 
                hist::Distribution2D&& p_wx, 
                hist::Distribution1D&& p_ww,
                hist::Distribution1D&& p_tot,
                double avg_displaced_V
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
            CompositeDistanceHistogramCrysol(
                hist::Distribution3D&& p_aa, 
                hist::Distribution3D&& p_ax, 
                hist::Distribution3D&& p_xx, 
                hist::Distribution2D&& p_aw, 
                hist::Distribution2D&& p_wx, 
                hist::Distribution1D&& p_ww, 
                hist::WeightedDistribution1D&& p_tot,
                double avg_displaced_V
            );

            const form_factor::storage::atomic::table_t& get_ff_table() const override;
            const form_factor::storage::cross::table_t& get_ffax_table() const override;
            const form_factor::storage::exv::table_t& get_ffxx_table() const override;

            Limit get_excluded_volume_scaling_factor_limits() const override;

            double average_displaced_V = 0;

        protected:
            double exv_factor(double q) const override;
            void initialize();

            inline static form_factor::storage::atomic::table_t ffaa_table;
            inline static form_factor::storage::cross::table_t  ffax_table;
            inline static form_factor::storage::exv::table_t    ffxx_table;
    };
}