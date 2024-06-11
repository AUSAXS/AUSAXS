#pragma once

#include "dataset/SimpleDataset.h"
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
     * @brief A class containing a single total distance histogram for use in the Debye transform.
     *        This is the simplest implementation of the transform; for more advanced features, see its subclasses.
     */
    class DistanceHistogram : protected Histogram {
        public: 
            /**
             * @brief Default constructor.
             */
            DistanceHistogram();

            /**
             * @brief Move constructor.
             */
            DistanceHistogram(DistanceHistogram&& other);

            /**
             * @brief Create an unweighted distance histogram.
             */
            DistanceHistogram(hist::Distribution1D&& p_tot);

            /**
             * @brief Create a weighted distance histogram. 
             *        A custom sinc(x) lookup table based on the weights in @a p_tot will be calculated and used in the Debye transform.
             */
            DistanceHistogram(hist::WeightedDistribution1D&& p_tot);

            /**
             * @brief Create a distance histogram from a composite distance histogram.
             *        Only the total distance histogram will be extracted from the composite histogram.
             */
            DistanceHistogram(std::unique_ptr<ICompositeDistanceHistogram> cdh);

            virtual ~DistanceHistogram() override;

            /**
             * @brief Perform the Fourier transform through the Debye equation.
             */
            virtual ScatteringProfile debye_transform() const;

            /**
             * @brief Perform the Fourier transform through the Debye equation.
             *        If the given q-axis is within the range of the default q-axis, they will be interpolated for better efficiency. 
             *        If not, a size q*d sinc(x) lookup table will be calculated for every call to this function.
             *
             * @param q The q values at which to evaluate the scattering. 
             */
            virtual SimpleDataset debye_transform(const std::vector<double>& q) const;

            /**
             * @brief Get the distance axis describing the current histogram.
             */
            const std::vector<double>& get_d_axis() const;

            /**
             * @brief Get the q axis used in the Fourier transform. 
             */
            static const std::vector<double>& get_q_axis();

            /**
             * @brief Get the total histogram counts. Equivalent to get_counts().
             */
            virtual const std::vector<double>& get_total_counts() const;
            std::vector<double>& get_total_counts(); // @copydoc get_total_counts() const

            /**
             * @brief Determine if the structure represented by this histogram is highly ordered.
             *        This is done by estimating the number of peaks in the histogram and comparing it to the number of bins.             
             */
            bool is_highly_ordered() const;

        protected:
            std::vector<double> d_axis;                             // The distance axis.
            std::unique_ptr<table::DebyeTable> weighted_sinc_table; // The weighted sinc table
            bool use_weighted_table = false;                        // Whether to use the weighted sinc table

            /**
             * @brief Get the sinc(x) lookup table for the Debye transform.
             */
            observer_ptr<const table::DebyeTable> get_sinc_table() const;

            /**
             * @brief Use a weighted sinc table for the Debye transform.
             *        This defines the weighted_sinc_table member based on the current d_axis and sets use_weighted_table to true.
             */
            void use_weighted_sinc_table();

        private:
            void initialize();
            void initialize(std::vector<double>&& d_axis);
    };
}