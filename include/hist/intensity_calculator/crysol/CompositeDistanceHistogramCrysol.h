#pragma once

#include <hist/intensity_calculator/CompositeDistanceHistogramFFExplicit.h>
#include <hist/distribution/Distribution3D.h>
#include <hist/distribution/Distribution2D.h>
#include <hist/distribution/Distribution1D.h>
#include <hist/distribution/WeightedDistribution3D.h>
#include <hist/distribution/WeightedDistribution2D.h>
#include <hist/distribution/WeightedDistribution1D.h>

namespace hist {
    /**
     * @brief An alternative to CompositeDistanceHistogramFFExplicit that mimics the CRYSOL excluded volume fitting. 
     */
    class CompositeDistanceHistogramCrysol : public CompositeDistanceHistogramFFExplicit {
        public:
            using CompositeDistanceHistogramFFExplicit::CompositeDistanceHistogramFFExplicit;
            ~CompositeDistanceHistogramCrysol() override = default;

            ScatteringProfile debye_transform() const override; // @copydoc DistanceHistogram::debye_transform() const
            virtual SimpleDataset debye_transform(const std::vector<double>& q) const override; // @copydoc DistanceHistogram::debye_transform(const std::vector<double>&) const
            virtual ScatteringProfile get_profile_aa() const override; // @copydoc ICompositeDistanceHistogram::get_profile_aa() const
            virtual ScatteringProfile get_profile_aw() const override; // @copydoc ICompositeDistanceHistogram::get_profile_aw() const
            virtual ScatteringProfile get_profile_ww() const override; // @copydoc ICompositeDistanceHistogram::get_profile_ww() const
            virtual ScatteringProfile get_profile_xx() const override; // @copydoc ICompositeDistanceHistogramExv::get_profile_xx() const
            virtual ScatteringProfile get_profile_ax() const override; // @copydoc ICompositeDistanceHistogramExv::get_profile_ax() const
            virtual ScatteringProfile get_profile_wx() const override; // @copydoc ICompositeDistanceHistogramExv::get_profile_wx() const
            Limit get_excluded_volume_scaling_factor_limits() const override; // @copydoc ICompositeDistanceHistogramExv::get_excluded_volume_scaling_factor_limits() const

        protected:
            double exv_factor(double q) const override;
    };
}