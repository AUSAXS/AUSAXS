// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/distribution/DistributionFwd.h>

#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <hist/distribution/Distribution1D.h>
#include <hist/distribution/Distribution2D.h>
#include <hist/distribution/Distribution3D.h>
#include <constants/Constants.h>

#include <vector>

namespace ausaxs::hist {
    /**
     * @brief A class containing partial distance histograms for the different types of interactions and atomic types. 
     *        Beyond the functionality of CompositeDistanceHistogram, this class also uses individual form factors for each atomic type.
     *        An average form factor is used for the excluded volume.
     *        For more information, see CompositeDistanceHistogram.
     * 
     * @tparam FormFactorTableType The form factor lookup table to use.
     */
    template<typename FormFactorTableType>
    class CompositeDistanceHistogramFFAvgBase : public ICompositeDistanceHistogramExv {
        public: 
            CompositeDistanceHistogramFFAvgBase();
            CompositeDistanceHistogramFFAvgBase(const CompositeDistanceHistogramFFAvgBase&);
            CompositeDistanceHistogramFFAvgBase(CompositeDistanceHistogramFFAvgBase&&) noexcept;
            CompositeDistanceHistogramFFAvgBase& operator=(const CompositeDistanceHistogramFFAvgBase&);
            CompositeDistanceHistogramFFAvgBase& operator=(CompositeDistanceHistogramFFAvgBase&&) noexcept;

            /**
             * @brief Create an unweighted composite distance histogram with form factors.
             * 
             * @param p_aa The partial distance histogram for atom-atom interactions.
             * @param p_aw The partial distance histogram for atom-water interactions.
             * @param p_ww The partial distance histogram for water-water interactions.
             * @param p_tot The total distance histogram. This is only used for determining the maximum distance.
             */
            CompositeDistanceHistogramFFAvgBase(
                hist::Distribution3D&& p_aa, 
                hist::Distribution2D&& p_aw, 
                hist::Distribution1D&& p_ww,
                hist::Distribution1D&& p_tot
            );

            /**
             * @brief Create a weighted composite distance histogram with form factors.
             * 
             * @param p_aa The partial distance histogram for atom-atom interactions.
             * @param p_aw The partial distance histogram for atom-water interactions.
             * @param p_ww The partial distance histogram for water-water interactions.
             * @param p_tot The total distance histogram. This is only used to extract the bin centers. 
             */
            CompositeDistanceHistogramFFAvgBase(
                hist::Distribution3D&& p_aa, 
                hist::Distribution2D&& p_aw, 
                hist::Distribution1D&& p_ww, 
                hist::WeightedDistribution1D&& p_tot
            );

            virtual ~CompositeDistanceHistogramFFAvgBase() override;

            virtual ScatteringProfile debye_transform() const override;
            virtual SimpleDataset debye_transform(const std::vector<double>& q) const override;

            void apply_water_scaling_factor(double k) override;
            void apply_excluded_volume_scaling_factor(double k) override;
            void apply_solvent_density_scaling_factor(double k) override;
            void apply_atomic_debye_waller_factor(double B) override;
            void apply_exv_debye_waller_factor(double B) override;

            /**
             * @brief Get the partial distance histogram for atom-atom interactions.
             */
            const Distribution1D& get_aa_counts() const override;
            Distribution1D& get_aa_counts() override; // @copydoc get_aa_counts() const

            /**
             * @brief Get the partial distance histogram for atom-water interactions.
             */
            const Distribution1D& get_aw_counts() const override;
            Distribution1D& get_aw_counts() override; // @copydoc get_aw_counts() const

            /**
             * @brief Get the partial distance histogram for water-water interactions.
             */
            const Distribution1D& get_ww_counts() const override;
            Distribution1D& get_ww_counts() override; // @copydoc get_ww_counts() const

            /**
             * @brief Get the partial distance histogram for atom-atom interactions.
             */
            const Distribution3D& get_aa_counts_ff() const;
            Distribution3D& get_aa_counts_ff(); // @copydoc get_aa_counts_ff() const

            /**
             * @brief Get the partial distance histogram for atom-water interactions.
             */
            const Distribution2D& get_aw_counts_ff() const;
            Distribution2D& get_aw_counts_ff(); // @copydoc get_aw_counts_ff() const

            /**
             * @brief Get the partial distance histogram for water-water interactions.
             */
            const Distribution1D& get_ww_counts_ff() const;
            Distribution1D& get_ww_counts_ff(); // @copydoc get_ww_counts_ff() const

            const std::vector<double>& get_counts() const override;

            virtual ScatteringProfile get_profile_aa() const override;
            virtual ScatteringProfile get_profile_aw() const override;
            virtual ScatteringProfile get_profile_ww() const override;
            virtual ScatteringProfile get_profile_ax() const override;
            virtual ScatteringProfile get_profile_xx() const override;
            virtual ScatteringProfile get_profile_wx() const override;

            virtual const FormFactorTableType& get_ff_table() const = 0;

            /**
             * @brief Get the atomic Debye Waller factor for a given q and sigma value.
             */
            static double get_atomic_debye_waller_factor(double q, double sigma);

            /**
             * @brief Get the excluded volume Debye Waller factor for a given q and sigma value.
             */
            static double get_exv_debye_waller_factor(double q, double sigma);

        protected:
            struct {
                double cw = 1;               // water density scaling factor
                double cx = 1;               // excluded volume scaling factor, method-dependent
                double crho = 1;             // solvent density scaling factor
                double DW_sigma_atomic = 0;  // atomic form factor debye-waller factor, zero for disabled
                double DW_sigma_exv = 0;     // excluded volume form factor debye-waller factor, zero for disabled
            } free_params;
            struct {Distribution3D aa; Distribution2D aw; Distribution1D ww;} distance_profiles;

            /**
             * @brief Get the q-dependent multiplicative factor for the excluded volume form factor.
             */
            virtual double exv_factor(double q) const;

        private:
            /**
             * @brief Get the atomic Debye Waller factor for a given q value.
             */
            double get_atomic_debye_waller_factor(double q) const;

            /**
             * @brief Get the excluded volume Debye Waller factor for a given q value.
             */
            double get_exv_debye_waller_factor(double q) const;

        //#################################//
        //###           CACHE           ###//
        //#################################//
        public:
            /**
             * @brief Get the cached intensity profiles.
             *        This may trigger a refresh if the cache is invalid.
             * 
             * @return [aa, ax, aw, xx, wx, ww]
             */
            [[nodiscard]] virtual std::tuple<
                std::vector<double>, std::vector<double>, std::vector<double>,
                std::vector<double>, std::vector<double>, std::vector<double> 
            > cache_get_intensity_profiles() const;

        protected:
            /**
             * @brief Get the cached total distance profiles. 
             *        This may trigger a refresh if the cache is invalid.
             * 
             * @return [aa, aw, ww]
             */
            [[nodiscard]] virtual std::tuple<const Distribution1D&, const Distribution1D&, const Distribution1D&> 
            cache_get_distance_profiles() const;

            mutable struct {
                // cached sinqd vals for each form factor combination
                // indexing as [ff1][ff2]
                mutable struct {
                    container::Container3D<double> aa;
                    container::Container2D<double> ax, aw;
                    container::Container1D<double> xx, wx, ww;
                    bool valid = false;
                } sinqd;

                mutable struct {
                    Distribution1D p_aa, p_aw, p_ww;
                    bool valid = false;
                } distance_profiles;

                mutable struct {
                    std::vector<double> aa, ax, aw, xx, wx, ww;
                    double cached_cx = -1;
                    double cached_cw = -1;
                    double cached_crho = -1;
                } intensity_profiles;
            } cache;

        private:
            /**
             * @brief Apply the Debye-Waller factors to the intensity profiles.
             */
            virtual std::tuple<
                std::vector<double>, std::vector<double>, std::vector<double>,
                std::vector<double>, std::vector<double>, std::vector<double> 
            > apply_debye_waller_factors(std::tuple<
                const std::vector<double>&, const std::vector<double>&, const std::vector<double>&,
                const std::vector<double>&, const std::vector<double>&, const std::vector<double>& 
            >) const;
            virtual void cache_refresh_intensity_profiles(bool sinqd_changed, bool cw_changed, bool cx_changed) const;
            virtual void cache_refresh_distance_profiles() const;
            virtual void cache_refresh_sinqd() const;
    };
}