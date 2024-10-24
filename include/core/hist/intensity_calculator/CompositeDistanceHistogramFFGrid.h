#pragma once

#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvg.h>
#include <form_factor/PrecalculatedFormFactorProduct.h>
#include <table/VectorDebyeTable.h>
#include <utility/TypeTraits.h>

namespace ausaxs::hist {
    /**
     * @brief A class containing partial distance histograms for the different types of interactions and atomic types. 
     *        Beyond the functionality of CompositeDistanceHistogram, this class also uses individual form factors for each atomic type.
     *        The excluded volume is approximated using a space-filling grid of spheres, filling the volume of the molecule. 
     *        This approach adds a substantial overhead to the calculations, but should give a more accurate representation of the excluded volume.
     */
    class CompositeDistanceHistogramFFGrid : public CompositeDistanceHistogramFFAvg {
        friend class CompositeDistanceHistogramFFGridSurface;
        friend class CompositeDistanceHistogramFFGridScalableExv;
        public:
            CompositeDistanceHistogramFFGrid(CompositeDistanceHistogramFFGrid&&) noexcept;
            CompositeDistanceHistogramFFGrid& operator=(CompositeDistanceHistogramFFGrid&&) noexcept;
            ~CompositeDistanceHistogramFFGrid() override;

            /**
             * @brief Create a weighted grid-based composite distance histogram with form factors.
             * 
             * @param p_aa The partial distance histogram for atom-atom interactions.
             * @param p_aw The partial distance histogram for atom-water interactions.
             * @param p_ww The partial distance histogram for water-water interactions.
             * @param p_tot_aa The total distance histogram for everything except the grid. This is only used to extract the bin centers.
             * @param p_tot_ax The total distance histogram for the cross terms. This is only used to extract the bin centers. 
             * @param p_tot_xx The total distance histogram for the grid only. Calculations involving this grid must use unique bin centers due to the highly ordered grid structure. 
             */
            CompositeDistanceHistogramFFGrid(
                hist::Distribution3D&& p_aa, 
                hist::Distribution2D&& p_aw, 
                hist::Distribution1D&& p_ww, 
                hist::WeightedDistribution1D&& p_tot_aa,
                hist::WeightedDistribution1D&& p_tot_ax,
                hist::WeightedDistribution1D&& p_tot_xx
            );

            const form_factor::storage::atomic::table_t& get_ff_table() const override {return ff_table;}

            /**
             * @brief Generate a new interal form factor table for the grid-based calculations.
             * 
             * @param ffx The excluded volume form factor to use. Leave as default to couple it to the grid volume.
             */
            template<FormFactorType T>
            static form_factor::storage::atomic::table_t generate_ff_table(T&& ffx);
            static form_factor::storage::atomic::table_t generate_ff_table(); //< @copydoc generate_ff_table(T&&)

            /**
             * @brief Regenerate the form factor table. This must be called to reflect changes in settings::grid::exv::width.
             * 
             * @param ffx The excluded volume form factor to use. Leave as default to couple it to the grid volume.
             */
            template<FormFactorType T>
            static void regenerate_ff_table(T&& ffx = {0});

            /**
             * @brief Get the distance axis for the excluded volume calculations. 
             *        If weighted bins are used, this will be distinct from the regular distance axis.
             */
            const std::vector<double>& get_d_axis_xx() const {return distance_axes.xx;}

            /**
             * @brief Get the distance axis for the cross term calculations. 
             *        If weighted bins are used, this will be distinct from the regular distance axis.
             */
            const std::vector<double>& get_d_axis_ax() const {return distance_axes.ax;}

        protected:
            /**
             * @brief Get the sinc(x) lookup table for the excluded volume for the Debye transform.
             */
            observer_ptr<const table::DebyeTable> get_sinc_table_xx() const;

            /**
             * @brief Get the sinc(x) lookup table for the cross terms for the Debye transform.
             */
            observer_ptr<const table::DebyeTable> get_sinc_table_ax() const;

            double exv_factor(double q) const override;

        private: 
            static form_factor::storage::atomic::table_t ff_table;
            struct {std::unique_ptr<table::VectorDebyeTable> xx, ax;} sinc_tables;
            struct {std::vector<double> xx, ax;} distance_axes;

            void initialize(std::vector<double>&& d_axis_ax, std::vector<double>&& d_axis_xx);

            //#################################//
            //###           CACHE           ###//
            //#################################//
            void cache_refresh_sinqd() const override;
    };
    static_assert(supports_nothrow_move_v<CompositeDistanceHistogramFFGrid>, "CompositeDistanceHistogramFFGrid should support nothrow move semantics.");
}