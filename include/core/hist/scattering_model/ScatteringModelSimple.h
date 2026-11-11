// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/distribution/Distribution1D.h>
#include <hist/scattering_model/IScatteringModel.h>

namespace ausaxs::hist {
    /**
     * @brief A scattering model where all scatterers are equal and have the same Gaussian form factor. 
     */
    class ScatteringModelSimple : public IScatteringModel {
        public: 
            ScatteringModelSimple(const ScatteringModelSimple&) = delete;
            ScatteringModelSimple(ScatteringModelSimple&&) noexcept;
            ScatteringModelSimple& operator=(ScatteringModelSimple&&) noexcept;
            ~ScatteringModelSimple() override;

            /**
             * @brief Create an unweighted scattering model.
             * 
             * @param p_aa The partial distance histogram for atom-atom interactions.
             * @param p_aw The partial distance histogram for atom-water interactions.
             * @param p_ww The partial distance histogram for water-water interactions.
             * @param p_tot The total distance histogram.
             */
            ScatteringModelSimple(
                hist::Distribution1D&& p_aa, 
                hist::Distribution1D&& p_aw, 
                hist::Distribution1D&& p_ww, 
                hist::Distribution1D&& p_tot
            );

            ScatteringModelSimple(
                hist::Distribution1D&& p_aa, 
                hist::Distribution1D&& p_aw, 
                hist::Distribution1D&& p_ww, 
                hist::WeightedDistribution1D&& p_tot
            );

            const hist::Distribution1D& get_aa_counts() const override;
            hist::Distribution1D& get_aa_counts() override;
            const Distribution1D& get_aw_counts() const override;
            Distribution1D& get_aw_counts() override;
            const Distribution1D& get_ww_counts() const override;
            Distribution1D& get_ww_counts() override;

            void apply_water_scaling_factor(double k) override;
            ScatteringProfile debye_transform() const override;
            using IScatteringModel::debye_transform;
            ScatteringProfile get_profile_aa() const override;
            ScatteringProfile get_profile_aw() const override;
            ScatteringProfile get_profile_ww() const override;

        protected:
            std::vector<double> debye_transform_direct(const std::vector<double>& q) const override;

        private:
            mutable struct {Distribution1D aa, aw, ww;} distance_profiles;

            ScatteringProfile transform_counts_gaussian(std::span<const double> counts) const;
    };
    static_assert(supports_nothrow_move_v<ScatteringModelSimple>, "ScatteringModelSimple should support nothrow move semantics.");
}