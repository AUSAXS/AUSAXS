#pragma once

#include <utility/Axis.h>

#include <vector>
#include <memory>

namespace dataset {class SimpleDataset;}
namespace hist {
    /**
     * @brief A class containing a single distance histogram.
     */
    class CompositeDistanceHistogram;
    class DistanceHistogram {
        public: 
            DistanceHistogram() = default;

            DistanceHistogram(std::vector<double>&& p_tot, const Axis& axis);

            /**
             * @brief Extract the total histogram from a CompositeDistanceHistogram.
             */
            DistanceHistogram(CompositeDistanceHistogram&& cdh);

            ~DistanceHistogram() = default;

            virtual Histogram debye_transform() const;

            virtual Histogram debye_transform(const std::vector<double>& q) const;

            std::vector<double>& get_total_histogram();

            std::vector<double> get_axis_vector() const;

            Axis& get_axis();

        private:
            std::vector<double> p_tot;
            Axis axis;
    };

    /**
     * @brief A class containing multiple partial distance histograms.
     */
    class CompositeDistanceHistogram : public DistanceHistogram {
        public: 
            CompositeDistanceHistogram() = default;

            CompositeDistanceHistogram(std::vector<double>&& p_pp, std::vector<double>&& p_hh, std::vector<double>&& p_hp, std::vector<double>&& p_tot, const Axis& axis);

            ~CompositeDistanceHistogram() = default;

            Histogram debye_transform() const override;

            Histogram debye_transform(const std::vector<double>& q) const override;

            std::vector<double>& get_pp_histogram();

            std::vector<double>& get_hh_histogram();

            std::vector<double>& get_hp_histogram();            

            void apply_water_scaling_factor(double k);

        private:
            std::vector<double> p_pp;
            std::vector<double> p_hp;
            std::vector<double> p_hh;
    };

    /**
     * @brief A class containing multiple partial distance histograms for multiple form factors. 
     */
    class CompositeDistanceHistogramFF : public DistanceHistogram {
        public: 
            CompositeDistanceHistogramFF() = default;

            ~CompositeDistanceHistogramFF() = default;

            Histogram debye_transform() const override;

            Histogram debye_transform(const std::vector<double>& q) const override;

            void apply_water_scaling_factor(double k);

        private:
            std::vector<CompositeDistanceHistogram> cdhs;
    };
}