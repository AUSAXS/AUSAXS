#pragma once

#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvg.h>
#include <form_factor/PrecalculatedFormFactorProduct.h>
#include <table/VectorDebyeTable.h>

namespace hist {
    /**
     * @brief A class containing partial distance histograms for the different types of interactions and atomic types. 
     *        Beyond the functionality of CompositeDistanceHistogram, this class also uses individual form factors for each atomic type.
     *        The excluded volume is approximated using a space-filling grid of spheres, filling the volume of the molecule. 
     *        This approach adds a substantial overhead to the calculations, but should give a more accurate representation of the excluded volume.
     */
    class CompositeDistanceHistogramFFGrid : public CompositeDistanceHistogramFFAvg {
        public:
            /**
             * @brief Create a grid-based unweighted composite distance histogram with form factors.
             * 
             * @param p_aa The partial distance histogram for atom-atom interactions.
             * @param p_aw The partial distance histogram for atom-water interactions.
             * @param p_ww The partial distance histogram for water-water interactions.
             * @param p_tot The total distance histogram. This is only used for determining the maximum distance.
             */
            CompositeDistanceHistogramFFGrid(
                hist::Distribution3D&& p_aa, 
                hist::Distribution2D&& p_aw, 
                hist::Distribution1D&& p_ww,
                hist::Distribution1D&& p_tot
            );

            /**
             * @brief Create a weighted grid-based composite distance histogram with form factors.
             * 
             * @param p_aa The partial distance histogram for atom-atom interactions.
             * @param p_aw The partial distance histogram for atom-water interactions.
             * @param p_ww The partial distance histogram for water-water interactions.
             * @param p_tot The total distance histogram for everything except the grid. This is only used to extract the bin centers.
             * @param p_tot_x The total distance histogram for the grid only. Calculations involving this grid must use unique bin centers due to the highly ordered grid structure. 
             */
            CompositeDistanceHistogramFFGrid(
                hist::Distribution3D&& p_aa, 
                hist::Distribution2D&& p_aw, 
                hist::Distribution1D&& p_ww, 
                hist::WeightedDistribution1D&& p_tot,
                hist::WeightedDistribution1D&& p_tot_x
            );

            /**
             * @brief Get the form factor table for the grid-based calculations.
             */
            const form_factor::storage::atomic::table_t& get_ff_table() const override {return ff_table;}

            /**
             * @brief Regenerate the form factor table.
             *        This is only necessary if the excluded volume radius has changed.
             */
            static void regenerate_table() {ff_table = generate_table();}

            virtual ScatteringProfile debye_transform() const override; // @copydoc DistanceHistogram::debye_transform() const

            /**
             * @brief Get the intensity profile for atom-water interactions.
             */
            virtual ScatteringProfile get_profile_xx() const override;

        private: 
            static form_factor::storage::atomic::table_t generate_table();
            inline static form_factor::storage::atomic::table_t ff_table = generate_table();
            std::unique_ptr<table::VectorDebyeTable> weighted_sinc_table_x;

            double exv_factor(double q) const override;

            /**
             * @brief Get the sinc(x) lookup table for the excluded volume for the Debye transform.
             */
            observer_ptr<const table::DebyeTable> get_sinc_table_x() const;

            void initialize(std::vector<double>&& d_axis_x);
    };
}