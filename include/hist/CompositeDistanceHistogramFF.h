#pragma once

#include <hist/CompositeDistanceHistogram.h>
#include <utility/Container1D.h>
#include <utility/Container2D.h>
#include <utility/Container3D.h>

#include <vector>

class Protein;
namespace hist {
    /**
     * @brief A class containing multiple partial distance histograms for multiple form factors. 
     */
    class CompositeDistanceHistogramFF : public CompositeDistanceHistogram {
        public: 
            CompositeDistanceHistogramFF() = default;

            CompositeDistanceHistogramFF(CompositeDistanceHistogramFF&& other) noexcept;

            /**
             * @brief Construct a new Composite Distance Histogram FF object
             * 
             * @param p_pp Partial distance histogram for atom-atom interactions
             * @param p_hp Partial distance histogram for atom-water interactions
             * @param p_hh Partial distance histogram for water-water interactions
             * @param p_tot Total distance histogram
             * @param axis Distance axis
             * @param Z_exv_avg Average charge of excluded volume
             */
            CompositeDistanceHistogramFF(Container3D<double>&& p_pp, Container2D<double>&& p_hp, Container1D<double>&& p_hh, std::vector<double>&& p_tot, const Axis& axis);

            ~CompositeDistanceHistogramFF() override;

            ScatteringProfile debye_transform() const override;

            SimpleDataset debye_transform(const std::vector<double>& q) const override;

            void apply_water_scaling_factor(double k) override;

            void apply_excluded_volume_scaling_factor(double k);

            const std::vector<double>& get_pp_counts() const override;

            const std::vector<double>& get_hh_counts() const override;

            const std::vector<double>& get_hp_counts() const override;

            const std::vector<double>& get_counts() const override;

        protected:
            double w_scaling = 1;
            double exv_scaling = 1;

        private:
            Container3D<double> cp_pp;
            Container2D<double> cp_hp;
            Container1D<double> cp_hh;
    };
}