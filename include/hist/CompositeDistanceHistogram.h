#pragma once

#include <hist/DistanceHistogram.h>

#include <vector>

namespace hist {
    class Histogram;

    /**
     * @brief A class containing multiple partial distance histograms.
     */
    class CompositeDistanceHistogram : public DistanceHistogram {
        public: 
            CompositeDistanceHistogram() = default;

            CompositeDistanceHistogram(std::vector<double>&& p_pp, std::vector<double>&& p_hp, std::vector<double>&& p_hh, std::vector<double>&& p_tot, const Axis& axis);

            virtual ~CompositeDistanceHistogram() override;

            const std::vector<double>& get_pp_histogram() const;
            const std::vector<double>& get_hh_histogram() const;
            const std::vector<double>& get_hp_histogram() const;            

            virtual void apply_water_scaling_factor(double k);

        private:
            std::vector<double> p_pp;
            std::vector<double> p_hp;
            std::vector<double> p_hh;
    };
}