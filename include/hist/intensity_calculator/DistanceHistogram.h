#pragma once

#include <utility/Axis.h>
#include <hist/distribution/Distribution1D.h>
#include <hist/Histogram.h>
#include <hist/HistFwd.h>
#include <table/ArrayDebyeTable.h>
#include <dataset/DatasetFwd.h>
#include <constants/Constants.h>

#include <vector>
#include <memory>

namespace hist {
    /**
     * @brief A DistanceHistogram is just a (x, count(x)) histogram.
     */
    class DistanceHistogram : public Histogram {
        public: 
            DistanceHistogram();
            DistanceHistogram(hist::Distribution1D&& p_tot, const Axis& axis);
            DistanceHistogram(hist::WeightedDistribution1D&& p_tot, const Axis& axis);

            /**
             * @brief Extract the total histogram from a CompositeDistanceHistogram.
             */
            DistanceHistogram(std::unique_ptr<ICompositeDistanceHistogram> cdh);

            virtual ~DistanceHistogram() override;

            virtual ScatteringProfile debye_transform() const;

            // virtual SimpleDataset debye_transform(const std::vector<double>& q) const;

            const std::vector<double>& get_d_axis() const;

            const std::vector<double>& get_q_axis() const;

            /**
             * @brief Get the total histogram counts. Equivalent to get_counts().
             */
            virtual const std::vector<double>& get_total_counts() const;

            /**
             * @brief Get the total histogram counts. Equivalent to get_counts().
             */
            std::vector<double>& get_total_counts();

        protected:
            std::vector<double> d_axis;                 // the distance axis
            std::vector<double> q_axis;                 // the q axis

            /**
             * @brief Get the sinc(x) lookup table for the Debye transform.
             */
            const table::ArrayDebyeTable& get_sinc_table() const;

            /**
             * @brief Use a weighted sinc table for the Debye transform.
             *        The weights are extracted from the WeightedHistogram struct, which automatically keeps track of all WeightedDistribution counts. 
             */
            void use_weighted_sinc_table();

        private:
            table::ArrayDebyeTable weighted_sinc_table; // the weighted sinc table
            bool use_weighted_table = false;            // whether to use the weighted sinc table

            void initialize();
    };
}