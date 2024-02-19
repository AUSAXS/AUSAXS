#pragma once

#include "hist/distribution/WeightedDistribution1D.h"
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/distribution/GenericDistribution2D.h>
#include <hist/distribution/GenericDistribution3D.h>
#include <constants/Constants.h>

#include <vector>

namespace hist {
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
            /**
             * @brief Default constructor.
             */
            CompositeDistanceHistogramFFAvgBase();

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

            // @copydoc DistanceHistogram::debye_transform() const
            virtual ScatteringProfile debye_transform() const override;

            // @copydoc DistanceHistogram::debye_transform(const std::vector<double>&) const
            virtual SimpleDataset debye_transform(const std::vector<double>& q) const override;

            /**
             * @brief Apply a scaling factor to the water partial distance histogram.
             */
            void apply_water_scaling_factor(double k) override;

            /**
             * @brief Apply a scaling factor to the excluded volume partial distance histogram.
             */
            void apply_excluded_volume_scaling_factor(double k) override;

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

            /**
             * @brief Get the total distance histogram.
             */
            const std::vector<double>& get_counts() const override;

            /**
             * @brief Get the intensity profile for atom-atom interactions.
             */
            virtual ScatteringProfile get_profile_aa() const override;

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
            virtual ScatteringProfile get_profile_ax() const override;

            /**
             * @brief Get the intensity profile for atom-water interactions.
             */
            virtual ScatteringProfile get_profile_xx() const override;

            /**
             * @brief Get the intensity profile for water-water interactions.
             */
            virtual ScatteringProfile get_profile_wx() const override;

            virtual const FormFactorTableType& get_ff_table() const = 0;

        protected:
            double cw = 1; // water scaling factor
            double cx = 1; // excluded volume scaling factor
            Distribution3D cp_aa; 
            Distribution2D cp_aw; 
            Distribution1D cp_ww;

        private:
            mutable Distribution1D p_aa;
            mutable Distribution1D p_aw;
            mutable Distribution1D p_ww;

            /**
             * @brief Get the q-dependent multiplicative factor for the excluded volume form factor.
             */
            virtual double exv_factor(double q) const;
    };
}