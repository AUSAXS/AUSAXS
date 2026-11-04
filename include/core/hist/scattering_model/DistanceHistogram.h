// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <constants/Constants.h>
#include <dataset/DatasetFwd.h>
#include <hist/HistFwd.h>
#include <hist/Histogram.h>
#include <hist/ScatteringProfile.h>
#include <hist/distribution/DistributionFwd.h>
#include <table/DebyeTableManager.h>
#include <utility/Axis.h>
#include <utility/TypeTraits.h>
#include <utility/observer_ptr.h>

#include <memory>
#include <span>
#include <vector>

namespace ausaxs::hist {
    /**
     * @brief A class containing a single total distance histogram for use in the Debye transform.
     *        This is the simplest implementation of the transform; for more advanced features, see its subclasses.
     */
    class DistanceHistogram : protected Histogram {
        public: 
            DistanceHistogram();
            DistanceHistogram(const DistanceHistogram&) = delete;
            DistanceHistogram(DistanceHistogram&& other) noexcept;
            DistanceHistogram& operator=(const DistanceHistogram&) = delete;
            DistanceHistogram& operator=(DistanceHistogram&& other) noexcept;

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
             * @brief Create a distance histogram from a scattering model.
             *        Only the total distance histogram will be extracted from the scattering model.
             */
            DistanceHistogram(std::unique_ptr<IScatteringModel> model);

            ~DistanceHistogram() override;

            /**
             * @brief Perform the Fourier transform through the Debye equation.
             *        Every scatterer carries the Gaussian form factor exp(-q²/2), so the intensity is multiplied by exp(-q²).
             *        Subclasses with their own form factors override this.
             */
            virtual ScatteringProfile debye_transform() const;

            /**
             * @brief Perform the Fourier transform through the Debye equation.
             *        q values within the default q axis (constants::axes::q_axis) are interpolated from the profile of this histogram.
             *        If settings::axes::qmin..qmax does not cover them, the q range is first expanded to do so; see expand_q_range.
             *        Outside the default q axis, a size q*d sinc(x) lookup table is calculated for every call to this function.
             *
             * @param q The q values at which to evaluate the scattering.
             */
            Dataset debye_transform(const std::vector<double>& q) const;

            /**
             * @brief Perform the Fourier transform through the Debye equation, optionally without the form factor.
             *
             * @tparam form_factor If true, this is identical to the virtual debye_transform().
             *                     If false, the plain Debye sum of the total histogram is returned, i.e. the intensity of point scatterers.
             */
            template<bool form_factor>
            ScatteringProfile debye_transform() const;

            /**
             * @brief Perform the Fourier transform through the Debye equation at the given q values, optionally without the form factor.
             *
             * @tparam form_factor If true, this is identical to debye_transform(q).
             *                     If false, the plain Debye sum of the total histogram is returned, i.e. the intensity of point scatterers.
             * @param q The q values at which to evaluate the scattering.
             */
            template<bool form_factor>
            Dataset debye_transform(const std::vector<double>& q) const;

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
            virtual const std::vector<double>& get_weighted_counts() const;

            /**
             * @brief Get the radius of gyration of the scattering described by this histogram.
             *
             * This is the q -> 0 limit of the Guinier expansion of debye_transform(), ln I(q) = ln I(0) - q^2 Rg^2/3 + O(q^4),
             * i.e. what a Guinier analysis of a perfect measurement of this model would return. Every scatterer therefore
             * enters with its excess form factor and the currently applied free parameters: the excluded volume is
             * subtracted, and the hydration shell counts with its contrast (cw), not with its raw electron count.
             *
             * The limit is extrapolated from the four lowest q-values of the transform. It is exact to 1e-4 A at the
             * default qmin, and degrades as qmin*Rg grows.
             *
             * @return The radius of gyration in A, or NaN if it is undefined: no contrast (I <= 0), or a negative Rg^2.
             */
            double get_Rg() const;

            /**
             * @brief Determine if the structure represented by this histogram is highly ordered.
             *        This is done by estimating the number of peaks in the histogram and comparing it to the number of bins.             
             */
            bool is_highly_ordered() const;

            /**
             * @brief Determine whether a distance-count sequence is highly ordered.
             *
             * This overload is used for partial histograms, where the atom-atom
             * component must be inspected instead of the total histogram.
             */
            static bool is_highly_ordered(const std::vector<double>& counts);

        protected:
            /**
             * @brief The Debye sum of @a counts at the given q values, optionally multiplied by the exp(-q²) form factor.
             *
             * @param sinqd_table The sinc table. Its row @a first_row + i must hold the sinc values for q[i].
             */
            template<bool form_factor>
            static std::vector<double> debye_sum(
                std::span<const double> counts, observer_ptr<const table::DebyeTable> sinqd_table, std::span<const double> q, int first_row = 0
            );

            /**
             * @brief The Debye sum of @a counts on the default q axis, optionally multiplied by the exp(-q²) form factor.
             *
             * @param sinqd_table A sinc table over the full default q axis.
             */
            template<bool form_factor>
            static ScatteringProfile debye_sum(std::span<const double> counts, observer_ptr<const table::DebyeTable> sinqd_table);

            /**
             * @brief Expand settings::axes::qmin..qmax to cover @a q, which must lie within the default q axis.
             *        Subclasses caching anything sized by the q range must invalidate it.
             */
            virtual void expand_q_range(const std::vector<double>& q) const;

            std::vector<double> d_axis;                     // The distance axis.
            mutable table::DebyeTableManager sinc_table;    // The sinc(x) lookup table manager for the Debye transform.

        private:
            void initialize();
            void initialize(std::vector<double>&& d_axis);
    };
    static_assert(supports_nothrow_move_v<DistanceHistogram>, "DistanceHistogram should support nothrow move semantics.");
}