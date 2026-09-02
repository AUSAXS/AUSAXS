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
     * @brief Common base for the excluded volume intensity calculators.
     *
     * The scattering amplitude is decomposed into an atomic, an excluded volume, and a hydration contribution,
     *      A = A_a - A_x + A_w,
     * so the intensity splits into the six partial profiles cached by this class and summed in debye_transform():
     *      I = aa + xx + ww + ax + xa + wx + xw + aw + wa
     * Note that p_aa already stores each atom-atom pair in both orderings, so its factor is baked into the counts, whereas p_aw stores each water-atom 
     * pair only once and needs the factor written out explicitly.
     *
     * The excluded volume belongs to the molecule only: A_x never extends over the hydration shell, whether it is built from dummy atoms placed on the 
     * real atoms or from a separate grid. There is therefore no water-exv self term, 
     * and wx has the single form factor pairing f_w*fx_a rather than a second pairing with an excluded volume of the water itself. 
     * This is valid for two reasons:
     *  (i) The hydration density is rescaled regardless. Instead of modelling the excess density as \delta\rho_w = (\rho_w - \rho_b), we model it as 
     *      \delta\rho_w = c_w*\rho'_w with c_w fitted. Our hydration shell is not physically accurate and does not reproduce the number density of a 
     *      real shell in the first place, so folding the bulk subtraction into the fitted scaling costs nothing.
     *
     *  (ii) It is a significant optimization. If the excluded volume also sat on the water molecules, any fit of c_x would depend on the generated shell, 
     *       which is stochastic; evaluating it properly would then require averaging over many generated shells.
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
             * @brief Get the raw (unweighted) partial distance histogram for atom-atom interactions, indexed by form factor type.
             *        These are the absolute distance counts before any form factor weighting.
             */
            const Distribution3D& get_raw_aa_counts_by_ff() const override;
            Distribution3D& get_raw_aa_counts_by_ff() override;

            /**
             * @brief Get the raw (unweighted) partial distance histogram for atom-water interactions, indexed by form factor type.
             *        These are the absolute distance counts before any form factor weighting.
             */
            const Distribution2D& get_raw_aw_counts_by_ff() const override;
            Distribution2D& get_raw_aw_counts_by_ff() override;

            /**
             * @brief Get the raw (unweighted) partial distance histogram for water-water interactions.
             *        These are the absolute distance counts before any form factor weighting.
             */
            const Distribution1D& get_raw_ww_counts_by_ff() const override;
            Distribution1D& get_raw_ww_counts_by_ff() override;

            /**
             * @brief Get the weighted partial distance histogram for atom-atom interactions, indexed by form factor type.
             *        These counts are scaled by form factor products f(q=0)*f(q=0) for each ff combination.
             * @deprecated Use get_aa_counts_ff() instead for raw counts. This method exists for backwards compatibility.
             */
            const Distribution3D& get_aa_counts_by_ff() const;
            Distribution3D& get_aa_counts_by_ff(); // @copydoc get_aa_counts_ff() const

            /**
             * @brief Get the weighted partial distance histogram for atom-water interactions.
             */
            const Distribution2D& get_aw_counts_by_ff() const;
            Distribution2D& get_aw_counts_by_ff(); // @copydoc get_aw_counts_ff() const

            /**
             * @brief Get the weighted partial distance histogram for water-water interactions.
             */
            const Distribution1D& get_ww_counts_by_ff() const;
            Distribution1D& get_ww_counts_by_ff(); // @copydoc get_ww_counts_ff() const

            const std::vector<double>& get_weighted_counts() const override;
            const std::vector<double>& get_total_raw_counts() const override;

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
                // cached sinqd vals for each form factor combination indexing as [ff1][ff2].
                // only the atomic terms live here. Models placing their dummy atoms on top of the real atoms reuse
                // these for the excluded volume as well, whereas models where the excluded volume is a separate set
                // of scatterers own the additional distributions themselves
                mutable struct {
                    container::Container3D<double> aa;
                    container::Container2D<double> aw;
                    container::Container1D<double> ww;
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

            /**
             * @brief Submit the jobs filling the excluded volume half of the sinqd cache. Implementations should only submit to the global pool 
             *        and not wait on it; cache_refresh_sinqd waits once for both halves so that they may overlap. 
             */
            virtual void cache_refresh_sinqd_exv() const = 0;

            /**
             * @brief Add the excluded volume terms (ax, xx, wx) to the intensity profile cache.
             *        The corresponding profiles have already been zeroed by the caller.
             * @param cx The excluded volume scaling factor evaluated for every q value.
             */
            virtual void cache_refresh_intensity_exv(const std::vector<double>& cx, bool cw_changed, bool cx_changed) const = 0;

        private:
            // @copydoc cache_refresh_sinqd_exv
            void cache_refresh_sinqd_atomic() const;

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