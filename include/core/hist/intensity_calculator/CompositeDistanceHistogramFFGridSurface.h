#pragma once

#include <hist/intensity_calculator/CompositeDistanceHistogramFFGrid.h>
#include <form_factor/PrecalculatedFormFactorProduct.h>
#include <table/VectorDebyeTable.h>
#include <utility/TypeTraits.h>

namespace hist {
    /**
     * @brief A class containing partial distance histograms for the different types of interactions and atomic types. 
     *        Beyond the functionality of CompositeDistanceHistogramFFGrid, this class allows scaling the form factors of only the surface grid cells. 
     */
    class CompositeDistanceHistogramFFGridSurface : public CompositeDistanceHistogramFFAvg {
        public:
            struct XXContainer {
                XXContainer(unsigned int size) : interior(size), surface(size), cross(size) {}
                WeightedDistribution1D interior;
                WeightedDistribution1D surface;
                WeightedDistribution1D cross;
                XXContainer operator+=(const XXContainer& other);
            };

            struct AXContainer {
                AXContainer(unsigned int ff, unsigned int size) : interior(ff, size), surface(ff, size) {}
                WeightedDistribution2D interior;
                WeightedDistribution2D surface;
                AXContainer operator+=(const AXContainer& other);
            };

            struct WXContainer {
                WXContainer(unsigned int size) : interior(size), surface(size) {}
                WeightedDistribution1D interior;
                WeightedDistribution1D surface;
                WXContainer operator+=(const WXContainer& other);
            };

            CompositeDistanceHistogramFFGridSurface(CompositeDistanceHistogramFFGridSurface&&) noexcept;
            CompositeDistanceHistogramFFGridSurface& operator=(CompositeDistanceHistogramFFGridSurface&&) noexcept;
            ~CompositeDistanceHistogramFFGridSurface() override;

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
            CompositeDistanceHistogramFFGridSurface(
                hist::Distribution3D&& p_aa, 
                hist::Distribution2D&& p_aw, 
                hist::Distribution1D&& p_ww, 
                XXContainer&& xx,
                AXContainer&& ax,
                WXContainer&& wx,
                hist::WeightedDistribution1D&& p_tot_aa,
                hist::WeightedDistribution1D&& p_tot_ax,
                hist::WeightedDistribution1D&& p_tot_xx
            );

            const form_factor::storage::atomic::table_t& get_ff_table() const override {return ff_table;} // @copydoc CompositeDistanceHistogramFFAvgBase::get_ff_table() const

            static void regenerate_table(); // @copydoc CompositeDistanceHistogramFFGrid::regenerate_table()

            // @copydoc DistanceHistogram::debye_transform() const
            virtual ScatteringProfile debye_transform() const override;

            // @copydoc DistanceHistogram::debye_transform(const std::vector<double>&) const
            virtual SimpleDataset debye_transform(const std::vector<double>& q) const override;

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

            Limit get_excluded_volume_scaling_factor_limits() const override; // @copydoc ICompositeDistanceHistogramExv::get_excluded_volume_scaling_factor_limits() const

            /**
             * @brief Get the excluded volume scaling factor.
             *
             * @param cx The scaling factor for the excluded volume.
             * @param q The scattering vector.
             */
            static double exv_factor(double q, double cx);

        private: 
            static form_factor::storage::atomic::table_t ff_table;
            struct {std::unique_ptr<table::VectorDebyeTable> xx, ax;} sinc_tables;
            struct {std::vector<double> xx, ax;} distance_axes;

            struct {hist::Distribution1D xx_i, xx_s, xx_c, wx_i, wx_s; hist::Distribution2D ax_i, ax_s;} exv_distance_profiles;
            hist::Distribution1D evaluate_xx_profile(double cx) const;
            hist::Distribution1D evaluate_wx_profile(double cx) const;
            hist::Distribution2D evaluate_ax_profile(double cx) const;

            double exv_factor(double q) const override;

            /**
             * @brief Get the sinc(x) lookup table for the excluded volume for the Debye transform.
             */
            observer_ptr<const table::DebyeTable> get_sinc_table_xx() const;

            /**
             * @brief Get the sinc(x) lookup table for the cross terms for the Debye transform.
             */
            observer_ptr<const table::DebyeTable> get_sinc_table_ax() const;

            void initialize(std::vector<double>&& d_axis_ax, std::vector<double>&& d_axis_xx);

            //#################################//
            //###           CACHE           ###//
            //#################################//
            template<bool sinqd_changed, bool cw_changed, bool cx_changed>
            void cache_refresh_intensity_profiles() const;
            void cache_refresh_sinqd() const;
    };
    static_assert(supports_nothrow_move_v<CompositeDistanceHistogramFFGridSurface>, "CompositeDistanceHistogramFFGridSurface should support nothrow move semantics.");
}