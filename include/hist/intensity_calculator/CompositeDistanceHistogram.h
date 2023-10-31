#pragma once

#include <hist/intensity_calculator/DistanceHistogram.h>

#include <vector>

namespace hist {
    /**
     * @brief A class containing multiple partial distance histograms.
     */
    class CompositeDistanceHistogram : public DistanceHistogram {
        public: 
            CompositeDistanceHistogram() = default;

            CompositeDistanceHistogram(std::vector<double>&& p_tot, const Axis& axis);

            CompositeDistanceHistogram(std::vector<double>&& p_aa, std::vector<double>&& p_wa, std::vector<double>&& p_ww, std::vector<double>&& p_tot, const Axis& axis);

            virtual ~CompositeDistanceHistogram() override;

            /**
             * @brief Get the partial distance histogram for atom-atom interactions.
             */
            virtual const std::vector<double>& get_aa_counts() const;
            virtual std::vector<double>& get_aa_counts(); // @copydoc get_aa_counts() const

            /**
             * @brief Get the partial distance histogram for atom-water interactions.
             */
            virtual const std::vector<double>& get_aw_counts() const;
            virtual std::vector<double>& get_aw_counts(); // @copydoc get_aw_counts() const

            /**
             * @brief Get the partial distance histogram for water-water interactions.
             */
            virtual const std::vector<double>& get_ww_counts() const;
            virtual std::vector<double>& get_ww_counts(); // @copydoc get_ww_counts() const

            /**
             * @brief Apply a scaling factor to the water partial distance histogram.
             */
            virtual void apply_water_scaling_factor(double k);

            /**
             * @brief Reset the water scaling factor to 1.
             */
            void reset_water_scaling_factor();

            /**
             * @brief Get the intensity profile for atom-atom interactions.
             */
            virtual const ScatteringProfile get_profile_aa() const;

            /**
             * @brief Get the intensity profile for atom-water interactions.
             */
            virtual const ScatteringProfile get_profile_aw() const;

            /**
             * @brief Get the intensity profile for water-water interactions.
             */
            virtual const ScatteringProfile get_profile_ww() const;

        protected:
            mutable std::vector<double> p_aa;
            mutable std::vector<double> p_aw;
            mutable std::vector<double> p_ww;
    };
}