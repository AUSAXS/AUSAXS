#pragma once

#include <hist/DistanceHistogram.h>

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

            virtual const std::vector<double>& get_pp_counts() const;

            virtual const std::vector<double>& get_hh_counts() const;

            virtual const std::vector<double>& get_hp_counts() const;            

            virtual void apply_water_scaling_factor(double k);

            void reset_water_scaling_factor();

        protected:
            mutable std::vector<double> p_aa;
            mutable std::vector<double> p_wa;
            mutable std::vector<double> p_ww;
    };
}