#pragma once

#include <utility/Axis.h>
#include <hist/distribution/Distribution1D.h>
#include <hist/distribution/WeightedDistribution1D.h>
#include <hist/Histogram.h>
#include <hist/HistFwd.h>
#include <table/DebyeTable.h>
#include <dataset/DatasetFwd.h>
#include <constants/Constants.h>
#include <utility/observer_ptr.h>

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

            static const std::vector<double>& get_q_axis();

            /**
             * @brief Get the total histogram counts. Equivalent to get_counts().
             */
            virtual const std::vector<double>& get_total_counts() const;

            /**
             * @brief Get the total histogram counts. Equivalent to get_counts().
             */
            std::vector<double>& get_total_counts();

        protected:
            std::vector<double> d_axis; // the distance axis

            /**
             * @brief Get the sinc(x) lookup table for the Debye transform.
             */
            observer_ptr<const table::DebyeTable> get_sinc_table() const;

            /**
             * @brief Use a weighted sinc table for the Debye transform.
             *        The weights are extracted from the WeightedHistogram struct, which automatically keeps track of all WeightedDistribution counts. 
             */
            void use_weighted_sinc_table(const std::vector<double>& weights);

        private:
            std::unique_ptr<table::DebyeTable> weighted_sinc_table;     // the weighted sinc table
            bool use_weighted_table = false;                            // whether to use the weighted sinc table

            void initialize();
    };
}