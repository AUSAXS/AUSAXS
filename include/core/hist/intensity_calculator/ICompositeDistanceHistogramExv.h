// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <container/Container1D.h>
#include <container/Container2D.h>
#include <container/Container3D.h>

namespace ausaxs::hist {
    class ICompositeDistanceHistogramExv : public ICompositeDistanceHistogram {
        public:
            using ICompositeDistanceHistogram::ICompositeDistanceHistogram;
            ICompositeDistanceHistogramExv() = default;
            ICompositeDistanceHistogramExv(const ICompositeDistanceHistogramExv&) = default;
            ICompositeDistanceHistogramExv(ICompositeDistanceHistogramExv&&) noexcept = default;
            ICompositeDistanceHistogramExv& operator=(const ICompositeDistanceHistogramExv&) = default;
            ICompositeDistanceHistogramExv& operator=(ICompositeDistanceHistogramExv&&) noexcept = default;
            virtual ~ICompositeDistanceHistogramExv() = default;

            /**
             * @brief Apply a scaling factor to the excluded volume partial distance histogram.
             */
            virtual void apply_excluded_volume_scaling_factor(double k) = 0;

            /**
             * @brief Apply a scaling factor to the solvent density used for excluded volume calculations. 
             */
            virtual void apply_solvent_density_scaling_factor(double k) = 0;

            /**
             * @brief Apply a Debye Waller B-factor to the atomic form factors.
             */
            virtual void apply_atomic_debye_waller_factor(double B) = 0;

            /**
             * @brief Apply a Debye Waller B-factor to the excluded volume form factors.
             */
            virtual void apply_exv_debye_waller_factor(double B) = 0;

            /**
             * @brief Get the partial intensity profile for atom-(excluded volume) interactions.
             */
            virtual ScatteringProfile get_profile_ax() const = 0;

            /**
             * @brief Get the partial intensity profile for (excluded volume)-(excluded volume) interactions.
             */
            virtual ScatteringProfile get_profile_xx() const = 0;

            /**
             * @brief Get the partial intensity profile for water-(excluded volume) interactions.
             */
            virtual ScatteringProfile get_profile_wx() const = 0;

            /**
             * @brief Get the limits for the excluded volume scaling factor parameter. 
             *        This is intended to be used by the fitter to set correct limits. 
             */
            virtual Limit get_excluded_volume_scaling_factor_limits() const;

            /**
             * @brief Get the limits for the solvent density scaling factor parameter. 
             *        This is intended to be used by the fitter to set correct limits. 
             */
            virtual Limit get_solvent_density_scaling_factor_limits() const;

            /**
             * @brief Get the limits for the Debye-Waller factor parameter. 
             *        This is intended to be used by the fitter to set correct limits. 
             */
            virtual Limit get_debye_waller_factor_limits() const;

            /**
             * @brief Get the raw unweighted total counts histogram (i.e., without form factor weighting).
             */
            virtual const std::vector<double>& get_total_raw_counts() const = 0;
            std::vector<double>& get_total_raw_counts(); // @copydoc get_total_raw_counts() const

            /**
             * @brief Get the total raw counts. This is equivalent to get_total_raw_counts().
             */
            const std::vector<double>& get_raw_counts() const;
            std::vector<double>& get_raw_counts(); // @copydoc get_raw_counts() const

            /**
             * @brief Get the raw (unweighted) partial distance histogram for atom-atom interactions, indexed by form factor type.
             *        These are the absolute distance counts before any form factor weighting.
             */
            virtual const Distribution3D& get_raw_aa_counts_by_ff() const = 0;
            virtual Distribution3D& get_raw_aa_counts_by_ff() = 0;

            /**
             * @brief Get the raw (unweighted) partial distance histogram for atom-water interactions, indexed by form factor type.
             *        These are the absolute distance counts before any form factor weighting.
             */
            virtual const Distribution2D& get_raw_aw_counts_by_ff() const = 0;
            virtual Distribution2D& get_raw_aw_counts_by_ff() = 0;

            /**
             * @brief Get the raw (unweighted) partial distance histogram for water-water interactions.
             *        These are the absolute distance counts before any form factor weighting.
             */
            virtual const Distribution1D& get_raw_ww_counts_by_ff() const = 0;
            virtual Distribution1D& get_raw_ww_counts_by_ff() = 0;
    };
}