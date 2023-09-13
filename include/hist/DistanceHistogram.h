#pragma once

#include <utility/Axis.h>
#include <hist/Histogram.h>

#include <vector>
#include <memory>

namespace table {class DebyeLookupTable;}
namespace dataset {class SimpleDataset;}
namespace hist {
    class CompositeDistanceHistogram;

    /**
     * @brief A DistanceHistogram is just a (x, count(x)) histogram.
     */
    class DistanceHistogram : public Histogram {
        public: 
            DistanceHistogram() = default;

            DistanceHistogram(std::vector<double>&& p_tot, const Axis& axis);

            /**
             * @brief Extract the total histogram from a CompositeDistanceHistogram.
             */
            DistanceHistogram(CompositeDistanceHistogram&& cdh);

            virtual ~DistanceHistogram() override;

            virtual Histogram debye_transform() const;

            virtual Histogram debye_transform(const std::vector<double>& q) const;

            Vector<double>& get_total_histogram();

            const std::vector<double>& get_d_axis() const;

            const std::vector<double>& get_q_axis() const;

            const Axis& get_axis() const;

        private:
			std::unique_ptr<table::DebyeLookupTable> sinqd_table;   // Lookup-table for sin(qd)/qd values for the scattering histograms.
            std::vector<double> d_axis;                             // the distance axis
            std::vector<double> q_axis;                             // the q axis

            void initialize();
    };
}