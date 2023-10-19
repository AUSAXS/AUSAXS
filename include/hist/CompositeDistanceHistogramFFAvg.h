#pragma once

#include <hist/CompositeDistanceHistogram.h>
#include <container/Container1D.h>
#include <container/Container2D.h>
#include <container/Container3D.h>

#include <vector>

namespace hist {
    /**
     * @brief A class containing multiple partial distance histograms for multiple form factors. 
     */
    class CompositeDistanceHistogramFFAvg : public CompositeDistanceHistogram {
        public: 
            CompositeDistanceHistogramFFAvg();

            CompositeDistanceHistogramFFAvg(CompositeDistanceHistogramFFAvg&& other) noexcept;

            /**
             * @brief Construct a new Composite Distance Histogram FF object
             * 
             * @param p_aa Partial distance histogram for atom-atom interactions
             * @param p_aw Partial distance histogram for atom-water interactions
             * @param p_ww Partial distance histogram for water-water interactions
             * @param p_tot Total distance histogram
             * @param axis Distance axis
             */
            CompositeDistanceHistogramFFAvg(container::Container3D<double>&& p_aa, container::Container2D<double>&& p_wa, container::Container1D<double>&& p_ww, std::vector<double>&& p_tot, const Axis& axis);

            virtual ~CompositeDistanceHistogramFFAvg() override;

            virtual ScatteringProfile debye_transform() const override;

            // SimpleDataset debye_transform(const std::vector<double>& q) const override;

            void apply_water_scaling_factor(double k) override;

            void apply_excluded_volume_scaling_factor(double k);

            const std::vector<double>& get_pp_counts() const override;

            const std::vector<double>& get_hh_counts() const override;

            const std::vector<double>& get_hp_counts() const override;

            const std::vector<double>& get_counts() const override;

        protected:
            double cw = 1; // water scaling factor
            double cx = 1; // excluded volume scaling factor
            container::Container3D<double> cp_aa;
            container::Container2D<double> cp_wa;
            container::Container1D<double> cp_ww;
    };
}