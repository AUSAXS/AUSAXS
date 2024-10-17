#pragma once

#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <container/Container1D.h>
#include <container/Container2D.h>
#include <container/Container3D.h>

namespace hist {
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
             * @brief Get the intensity profile for atom-(excluded volume) interactions.
             */
            virtual ScatteringProfile get_profile_ax() const = 0;

            /**
             * @brief Get the intensity profile for (excluded volume)-(excluded volume) interactions.
             */
            virtual ScatteringProfile get_profile_xx() const = 0;

            /**
             * @brief Get the intensity profile for water-(excluded volume) interactions.
             */
            virtual ScatteringProfile get_profile_wx() const = 0;

            /**
             * @brief Get the limits for the excluded volume scaling factor parameter. 
             *        This is intended to be used by the fitter to set correct limits. 
             */
            virtual Limit get_excluded_volume_scaling_factor_limits() const {return {0.92, 1.08};}

            /**
             * @brief Get the limits for the solvent density scaling factor parameter. 
             *        This is intended to be used by the fitter to set correct limits. 
             */
            virtual Limit get_solvent_scaling_factor_limits() const {return {0.95, 1.05};}
    };
}