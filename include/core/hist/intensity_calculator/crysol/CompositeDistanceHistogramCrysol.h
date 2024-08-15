#pragma once

#include <hist/intensity_calculator/CompositeDistanceHistogramFFExplicit.h>
#include <hist/intensity_calculator/crysol/FormFactorCrysol.h>

namespace hist {
    /**
     * @brief An alternative to CompositeDistanceHistogramFFExplicit that mimics the CRYSOL excluded volume fitting. 
     */
    
    class CompositeDistanceHistogramCrysol : public CompositeDistanceHistogramFFExplicitBase<form_factor::storage::atomic::table_t, form_factor::storage::cross::table_t, form_factor::storage::exv::table_t>{
        public:
            using CompositeDistanceHistogramFFExplicitBase::CompositeDistanceHistogramFFExplicitBase;
            ~CompositeDistanceHistogramCrysol() override = default;

            const form_factor::storage::atomic::table_t& get_ff_table() const override;
            const form_factor::storage::cross::table_t& get_ffax_table() const override;
            const form_factor::storage::exv::table_t& get_ffxx_table() const override;

            Limit get_excluded_volume_scaling_factor_limits() const override;

            ScatteringProfile debye_transform() const override; // @copydoc DistanceHistogram::debye_transform() const
            virtual SimpleDataset debye_transform(const std::vector<double>& q) const override; // @copydoc DistanceHistogram::debye_transform(const std::vector<double>&) const
            virtual ScatteringProfile get_profile_aa() const override; // @copydoc ICompositeDistanceHistogram::get_profile_aa() const
            virtual ScatteringProfile get_profile_aw() const override; // @copydoc ICompositeDistanceHistogram::get_profile_aw() const
            virtual ScatteringProfile get_profile_ww() const override; // @copydoc ICompositeDistanceHistogram::get_profile_ww() const
            virtual ScatteringProfile get_profile_xx() const override; // @copydoc ICompositeDistanceHistogramExv::get_profile_xx() const
            virtual ScatteringProfile get_profile_ax() const override; // @copydoc ICompositeDistanceHistogramExv::get_profile_ax() const
            virtual ScatteringProfile get_profile_wx() const override; // @copydoc ICompositeDistanceHistogramExv::get_profile_wx() const

        protected:
            double exv_factor(double q) const override;
            inline static form_factor::storage::atomic::table_t ffaa_table = form_factor::crysol::storage::atomic::generate_table();
            inline static form_factor::storage::cross::table_t  ffax_table = form_factor::crysol::storage::cross::generate_table();
            inline static form_factor::storage::exv::table_t    ffxx_table = form_factor::crysol::storage::exv::generate_table();
    };
}