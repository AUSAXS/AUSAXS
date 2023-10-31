#pragma once

#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvg.h>
#include <container/Container1D.h>
#include <container/Container2D.h>
#include <container/Container3D.h>

#include <vector>

namespace hist {
    /**
     * @brief A class containing multiple partial distance histograms for multiple form factors. 
     */
    class CompositeDistanceHistogramFFExplicit : public CompositeDistanceHistogramFFAvg {
        public: 
            CompositeDistanceHistogramFFExplicit();

            /**
             * @brief Construct a new Composite Distance Histogram FF object
             * 
             * @param p_aa Partial distance histogram for atom-atom interactions
             * @param p_ax Partial distance histogram for atom-excluded volume interactions
             * @param p_xx Partial distance histogram for excluded volume-excluded volume interactions
             * @param p_aw Partial distance histogram for atom-water interactions
             * @param p_wx Partial distance histogram for water-excluded volume interactions
             * @param p_ww Partial distance histogram for water-water interactions
             * @param p_tot Total distance histogram
             * @param axis Distance axis
             */
            CompositeDistanceHistogramFFExplicit(
                container::Container3D<double>&& p_aa, container::Container3D<double>&& p_ax, container::Container3D<double>&& p_xx, 
                container::Container2D<double>&& p_wa, container::Container2D<double>&& p_wx, container::Container1D<double>&& p_ww, 
                std::vector<double>&& p_tot, const Axis& axis);

            ~CompositeDistanceHistogramFFExplicit() override;

            ScatteringProfile debye_transform() const override;

            /**
             * @brief Get the intensity profile for atom-atom interactions.
             */
            virtual const ScatteringProfile get_profile_ax() const;

            /**
             * @brief Get the intensity profile for atom-water interactions.
             */
            virtual const ScatteringProfile get_profile_xx() const;

            /**
             * @brief Get the intensity profile for water-water interactions.
             */
            virtual const ScatteringProfile get_profile_wx() const;

        protected:
            double G_factor(double q) const;

        private:
            container::Container3D<double> cp_ax;
            container::Container3D<double> cp_xx;
            container::Container2D<double> cp_wx;
    };
}