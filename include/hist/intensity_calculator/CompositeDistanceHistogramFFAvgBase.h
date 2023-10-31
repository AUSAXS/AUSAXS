#pragma once

#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <container/Container1D.h>
#include <container/Container2D.h>
#include <container/Container3D.h>

#include <vector>

namespace hist {
    /**
     * @brief A class containing multiple partial distance histograms for multiple form factors. 
     */
    template<typename T>
    class CompositeDistanceHistogramFFAvgBase : public CompositeDistanceHistogram {
        public: 
            CompositeDistanceHistogramFFAvgBase();

            /**
             * @brief Construct a new Composite Distance Histogram FF object
             * 
             * @param p_aa Partial distance histogram for atom-atom interactions
             * @param p_aw Partial distance histogram for atom-water interactions
             * @param p_ww Partial distance histogram for water-water interactions
             * @param p_tot Total distance histogram
             * @param axis Distance axis
             */
            CompositeDistanceHistogramFFAvgBase(container::Container3D<double>&& p_aa, container::Container2D<double>&& p_wa, container::Container1D<double>&& p_ww, std::vector<double>&& p_tot, const Axis& axis);

            virtual ~CompositeDistanceHistogramFFAvgBase() override;

            virtual ScatteringProfile debye_transform() const override;

            /**
             * @brief Apply a scaling factor to the water partial distance histogram.
             */
            void apply_water_scaling_factor(double k) override;

            /**
             * @brief Apply a scaling factor to the excluded volume partial distance histogram.
             */
            void apply_excluded_volume_scaling_factor(double k);

            /**
             * @brief Get the partial distance histogram for atom-atom interactions.
             */
            const std::vector<double>& get_aa_counts() const override;
            std::vector<double>& get_aa_counts() override; // @copydoc get_aa_counts() const

            /**
             * @brief Get the partial distance histogram for atom-water interactions.
             */
            const std::vector<double>& get_aw_counts() const override;
            std::vector<double>& get_aw_counts() override; // @copydoc get_aw_counts() const

            /**
             * @brief Get the partial distance histogram for water-water interactions.
             */
            const std::vector<double>& get_ww_counts() const override;
            std::vector<double>& get_ww_counts() override; // @copydoc get_ww_counts() const

            /**
             * @brief Get the partial distance histogram for atom-atom interactions.
             */
            const container::Container3D<double>& get_aa_counts_ff() const;
            container::Container3D<double>& get_aa_counts_ff(); // @copydoc get_aa_counts_ff() const

            /**
             * @brief Get the partial distance histogram for atom-water interactions.
             */
            const container::Container2D<double>& get_aw_counts_ff() const;
            container::Container2D<double>& get_aw_counts_ff(); // @copydoc get_aw_counts_ff() const

            /**
             * @brief Get the partial distance histogram for water-water interactions.
             */
            const container::Container1D<double>& get_ww_counts_ff() const;
            container::Container1D<double>& get_ww_counts_ff(); // @copydoc get_ww_counts_ff() const

            /**
             * @brief Get the total distance histogram.
             */
            const std::vector<double>& get_counts() const override;

            /**
             * @brief Get the intensity profile for atom-atom interactions.
             */
            virtual const ScatteringProfile get_profile_aa() const override;

            /**
             * @brief Get the intensity profile for atom-water interactions.
             */
            virtual const ScatteringProfile get_profile_aw() const override;

            /**
             * @brief Get the intensity profile for water-water interactions.
             */
            virtual const ScatteringProfile get_profile_ww() const override;

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

            virtual const T& get_ff_table() const = 0;

        protected:
            double cw = 1; // water scaling factor
            double cx = 1; // excluded volume scaling factor
            container::Container3D<double> cp_aa;
            container::Container2D<double> cp_aw;
            container::Container1D<double> cp_ww;
    };
}