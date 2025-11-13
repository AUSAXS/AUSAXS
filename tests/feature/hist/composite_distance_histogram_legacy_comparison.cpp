#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hist/histogram_manager/HistogramManagerMTFFAvg.h>
#include <hist/histogram_manager/HistogramManagerMTFFExplicit.h>
#include <hist/histogram_manager/HistogramManagerMTFFGrid.h>
#include <hist/histogram_manager/HistogramManagerMTFFGridSurface.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvg.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFExplicit.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGrid.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGridSurface.h>
#include <data/Molecule.h>
#include <grid/Grid.h>
#include <dataset/SimpleDataset.h>
#include <plots/All.h>
#include <settings/All.h>

#include "hist/hist_test_helper.h"

using namespace ausaxs;
using namespace ausaxs::hist;

/*
    This file tests the legacy single-threaded intensity calculations against the new cached multi-threaded implementations, ensuring they are consistent. 
*/

struct DebugCompositeDistanceHistogramFFAvg : public CompositeDistanceHistogramFFAvg {
    using CompositeDistanceHistogramFFAvg::CompositeDistanceHistogramFFAvg;
    using CompositeDistanceHistogramFFAvg::cache_get_intensity_profiles;
    DebugCompositeDistanceHistogramFFAvg(CompositeDistanceHistogramFFAvg&& other) : CompositeDistanceHistogramFFAvg(std::move(other)) {}

    inline static bool called_aa = false;
    inline static bool called_aw = false;
    inline static bool called_ww = false;
    inline static bool called_ax = false;
    inline static bool called_wx = false;
    inline static bool called_xx = false;
    inline static bool called_debye = false;
    ScatteringProfile get_profile_aa() const override  {
        called_aa = true;
        const auto& ff_table = get_ff_table();
        auto sinqd_table = sinc_table.get_sinc_table();
        Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
        unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

        std::vector<double> Iq(debye_axis.bins, 0);
        for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
            for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
                for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                    double aa_sum = std::inner_product(distance_profiles.aa.begin(ff1, ff2), distance_profiles.aa.end(ff1, ff2), sinqd_table->begin(q), 0.0);
                    Iq[q-q0] += aa_sum*ff_table.index(ff1, ff2).evaluate(q);
                }
            }
        }
        return ScatteringProfile(std::move(Iq), debye_axis);
    }

    ScatteringProfile get_profile_ax() const override  {
        called_ax = true;
        const auto& ff_table = get_ff_table();
        auto sinqd_table = sinc_table.get_sinc_table();
        Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
        unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

        std::vector<double> Iq(debye_axis.bins, 0);
        for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
            double cx = exv_factor(q);
            for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
                double ax_sum = std::inner_product(distance_profiles.aa.begin(ff1, form_factor::exv_bin), distance_profiles.aa.end(ff1, form_factor::exv_bin), sinqd_table->begin(q), 0.0);
                Iq[q-q0] += 2*cx*ax_sum*ff_table.index(ff1, form_factor::exv_bin).evaluate(q);
            }
        }
        return ScatteringProfile(std::move(Iq), debye_axis);
    }

    ScatteringProfile get_profile_xx() const override  {
        called_xx = true;
        const auto& ff_table = get_ff_table();
        auto sinqd_table = sinc_table.get_sinc_table();
        Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
        unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

        std::vector<double> Iq(debye_axis.bins, 0);
        for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
            double cx = exv_factor(q);
            double xx_sum = std::inner_product(distance_profiles.aa.begin(form_factor::exv_bin, form_factor::exv_bin), distance_profiles.aa.end(form_factor::exv_bin, form_factor::exv_bin), sinqd_table->begin(q), 0.0);
            Iq[q-q0] += cx*cx*xx_sum*ff_table.index(form_factor::exv_bin, form_factor::exv_bin).evaluate(q);
        }
        return ScatteringProfile(std::move(Iq), debye_axis);
    }

    ScatteringProfile get_profile_wx() const override  {
        called_wx = true;
        const auto& ff_table = get_ff_table();
        auto sinqd_table = sinc_table.get_sinc_table();
        Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
        unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

        std::vector<double> Iq(debye_axis.bins, 0);
        for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
            double cx = exv_factor(q);
            double ew_sum = std::inner_product(distance_profiles.aw.begin(form_factor::exv_bin), distance_profiles.aw.end(form_factor::exv_bin), sinqd_table->begin(q), 0.0);
            Iq[q-q0] += 2*cx*free_params.cw*ew_sum*ff_table.index(form_factor::exv_bin, form_factor::water_bin).evaluate(q);
        }
        return ScatteringProfile(std::move(Iq), debye_axis);
    }

    ScatteringProfile get_profile_aw() const override  {
        called_aw = true;
        const auto& ff_table = get_ff_table();
        auto sinqd_table = sinc_table.get_sinc_table();
        Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
        unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

        std::vector<double> Iq(debye_axis.bins, 0);
        for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
            for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
                double aw_sum = std::inner_product(distance_profiles.aw.begin(ff1), distance_profiles.aw.end(ff1), sinqd_table->begin(q), 0.0);
                Iq[q-q0] += 2*free_params.cw*aw_sum*ff_table.index(ff1, form_factor::water_bin).evaluate(q);
            }
        }
        return ScatteringProfile(std::move(Iq), debye_axis);
    }

    ScatteringProfile get_profile_ww() const override  {
        called_ww = true;
        const auto& ff_table = get_ff_table();
        auto sinqd_table = sinc_table.get_sinc_table();
        Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
        unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

        std::vector<double> Iq(debye_axis.bins, 0);
        for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
            double ww_sum = std::inner_product(distance_profiles.ww.begin(), distance_profiles.ww.end(), sinqd_table->begin(q), 0.0);
            Iq[q-q0] += free_params.cw*free_params.cw*ww_sum*ff_table.index(form_factor::water_bin, form_factor::water_bin).evaluate(q);
        }
        return ScatteringProfile(std::move(Iq), debye_axis);
    }

    ScatteringProfile debye_transform() const override {
        called_debye = true;
        const auto& ff_table = get_ff_table();
        auto sinqd_table = sinc_table.get_sinc_table();

        // calculate the Debye scattering intensity
        Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
        unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

        std::vector<double> Iq(debye_axis.bins, 0);
        for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
            double cx = exv_factor(q);
            for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
                // atom-atom
                for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                    double aa_sum = std::inner_product(distance_profiles.aa.begin(ff1, ff2), distance_profiles.aa.end(ff1, ff2), sinqd_table->begin(q), 0.0);
                    Iq[q-q0] += aa_sum*ff_table.index(ff1, ff2).evaluate(q);
                }

                // atom-exv
                double ax_sum = std::inner_product(distance_profiles.aa.begin(ff1, form_factor::exv_bin), distance_profiles.aa.end(ff1, form_factor::exv_bin), sinqd_table->begin(q), 0.0);
                Iq[q-q0] -= 2*cx*ax_sum*ff_table.index(ff1, form_factor::exv_bin).evaluate(q);

                // atom-water
                double aw_sum = std::inner_product(distance_profiles.aw.begin(ff1), distance_profiles.aw.end(ff1), sinqd_table->begin(q), 0.0);
                Iq[q-q0] += 2*free_params.cw*aw_sum*ff_table.index(ff1, form_factor::water_bin).evaluate(q);
            }

            // exv-exv
            double xx_sum = std::inner_product(distance_profiles.aa.begin(form_factor::exv_bin, form_factor::exv_bin), distance_profiles.aa.end(form_factor::exv_bin, form_factor::exv_bin), sinqd_table->begin(q), 0.0);
            Iq[q-q0] += cx*cx*xx_sum*ff_table.index(form_factor::exv_bin, form_factor::exv_bin).evaluate(q);

            // exv-water
            double ew_sum = std::inner_product(distance_profiles.aw.begin(form_factor::exv_bin), distance_profiles.aw.end(form_factor::exv_bin), sinqd_table->begin(q), 0.0);
            Iq[q-q0] -= 2*cx*free_params.cw*ew_sum*ff_table.index(form_factor::exv_bin, form_factor::water_bin).evaluate(q);

            // water-water
            double ww_sum = std::inner_product(distance_profiles.ww.begin(), distance_profiles.ww.end(), sinqd_table->begin(q), 0.0);
            Iq[q-q0] += free_params.cw*free_params.cw*ww_sum*ff_table.index(form_factor::water_bin, form_factor::water_bin).evaluate(q);
        }
        return ScatteringProfile(Iq, debye_axis);
    }
};

struct DebugCompositeDistanceHistogramFFExplicit : public CompositeDistanceHistogramFFExplicit {
    using CompositeDistanceHistogramFFExplicit::CompositeDistanceHistogramFFExplicit;
    using CompositeDistanceHistogramFFExplicit::cache_get_intensity_profiles;
    DebugCompositeDistanceHistogramFFExplicit(CompositeDistanceHistogramFFExplicit&& other) : CompositeDistanceHistogramFFExplicit(std::move(other)) {}

    inline static bool called_aa = false;
    inline static bool called_aw = false;
    inline static bool called_ww = false;
    inline static bool called_ax = false;
    inline static bool called_wx = false;
    inline static bool called_xx = false;
    inline static bool called_debye = false;
    ScatteringProfile get_profile_aa() const override  {
        called_aa = true;
        const auto& ff_table = get_ff_table();
        auto sinqd_table = sinc_table.get_sinc_table();
        Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
        unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

        std::vector<double> Iq(debye_axis.bins, 0);
        for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
            for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
                for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                    double aa_sum = std::inner_product(distance_profiles.aa.begin(ff1, ff2), distance_profiles.aa.end(ff1, ff2), sinqd_table->begin(q), 0.0);
                    Iq[q-q0] += aa_sum*ff_table.index(ff1, ff2).evaluate(q);
                }
            }
        }
        return ScatteringProfile(std::move(Iq), debye_axis);
    }

    ScatteringProfile get_profile_aw() const override  {
        called_aw = true;
        const auto& ff_table = get_ff_table();
        auto sinqd_table = sinc_table.get_sinc_table();
        Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
        unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

        std::vector<double> Iq(debye_axis.bins, 0);
        for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
            for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
                double aw_sum = std::inner_product(distance_profiles.aw.begin(ff1), distance_profiles.aw.end(ff1), sinqd_table->begin(q), 0.0);
                Iq[q-q0] += 2*free_params.cw*aw_sum*ff_table.index(ff1, form_factor::water_bin).evaluate(q);
            }
        }
        return ScatteringProfile(std::move(Iq), debye_axis);
    }

    ScatteringProfile get_profile_ww() const override  {
        called_ww = true;
        const auto& ff_table = get_ff_table();
        auto sinqd_table = sinc_table.get_sinc_table();
        Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
        unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

        std::vector<double> Iq(debye_axis.bins, 0);
        for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
            double ww_sum = std::inner_product(distance_profiles.ww.begin(), distance_profiles.ww.end(), sinqd_table->begin(q), 0.0);
            Iq[q-q0] += free_params.cw*free_params.cw*ww_sum*ff_table.index(form_factor::water_bin, form_factor::water_bin).evaluate(q);
        }
        return ScatteringProfile(std::move(Iq), debye_axis);
    }

    ScatteringProfile get_profile_ax() const override  {
        called_ax = true;
        const auto& ff_ax_table = get_ffax_table();
        auto sinqd_table = this->sinc_table.get_sinc_table();
        Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
        unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

        std::vector<double> Iq(debye_axis.bins, 0);
        for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
            double cx = exv_factor(constants::axes::q_vals[q]);
            for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
                for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                    double ax_sum = std::inner_product(exv_distance_profiles.ax.begin(ff1, ff2), exv_distance_profiles.ax.end(ff1, ff2), sinqd_table->begin(q), 0.0);
                    Iq[q-q0] += 2*cx*ax_sum*ff_ax_table.index(ff1, ff2).evaluate(q);
                }
            }
        }
        return ScatteringProfile(std::move(Iq), debye_axis);
    }

    ScatteringProfile get_profile_xx() const override  {
        called_xx = true;
        const auto& ff_xx_table = get_ffxx_table();
        auto sinqd_table = this->sinc_table.get_sinc_table();
        Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
        unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

        std::vector<double> Iq(debye_axis.bins, 0);
        for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
            double cx2 = std::pow(exv_factor(constants::axes::q_vals[q]), 2);
            for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
                for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                    double xx_sum = std::inner_product(exv_distance_profiles.xx.begin(ff1, ff2), exv_distance_profiles.xx.end(ff1, ff2), sinqd_table->begin(q), 0.0);
                    Iq[q-q0] += cx2*xx_sum*ff_xx_table.index(ff1, ff2).evaluate(q);
                }
            }
        }
        return ScatteringProfile(std::move(Iq), debye_axis);
    }

    ScatteringProfile get_profile_wx() const override  {
        called_wx = true;
        const auto& ff_ax_table = get_ffax_table();
        auto sinqd_table = this->sinc_table.get_sinc_table();
        Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
        unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

        std::vector<double> Iq(debye_axis.bins, 0);
        unsigned int ff_w_index = static_cast<int>(form_factor::form_factor_t::OH);
        for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
            double cx = exv_factor(constants::axes::q_vals[q]);
            for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
                double wx_sum = std::inner_product(exv_distance_profiles.wx.begin(ff1), exv_distance_profiles.wx.end(ff1), sinqd_table->begin(q), 0.0);
                Iq[q-q0] += 2*cx*this->free_params.cw*wx_sum*ff_ax_table.index(ff_w_index, ff1).evaluate(q);
            }
        }
        return ScatteringProfile(std::move(Iq), debye_axis);
    }

    ScatteringProfile debye_transform() const override  {
        called_debye = true;
        const auto& ff_aa_table = get_ffaa_table();
        const auto& ff_ax_table = get_ffax_table();
        const auto& ff_xx_table = get_ffxx_table();
        auto sinqd_table = this->sinc_table.get_sinc_table();

        // calculate the Debye scattering intensity
        Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
        unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

        std::vector<double> Iq(debye_axis.bins, 0);
        for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
            double cx = exv_factor(constants::axes::q_vals[q]);
            for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
                for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                    // atom-atom
                    double aa_sum = std::inner_product(this->distance_profiles.aa.begin(ff1, ff2), this->distance_profiles.aa.end(ff1, ff2), sinqd_table->begin(q), 0.0);
                    Iq[q-q0] += aa_sum*ff_aa_table.index(ff1, ff2).evaluate(q);

                    // atom-exv
                    double ax_sum = std::inner_product(exv_distance_profiles.ax.begin(ff1, ff2), exv_distance_profiles.ax.end(ff1, ff2), sinqd_table->begin(q), 0.0);
                    Iq[q-q0] -= 2*cx*ax_sum*ff_ax_table.index(ff1, ff2).evaluate(q);

                    // exv-exv
                    double xx_sum = std::inner_product(exv_distance_profiles.xx.begin(ff1, ff2), exv_distance_profiles.xx.end(ff1, ff2), sinqd_table->begin(q), 0.0);
                    Iq[q-q0] += cx*cx*xx_sum*ff_xx_table.index(ff1, ff2).evaluate(q);
                }

                // atom-water
                double aw_sum = std::inner_product(this->distance_profiles.aw.begin(ff1), this->distance_profiles.aw.end(ff1), sinqd_table->begin(q), 0.0);
                Iq[q-q0] += 2*this->free_params.cw*aw_sum*ff_aa_table.index(ff1, form_factor::water_bin).evaluate(q);

                // exv-water
                double wx_sum = std::inner_product(exv_distance_profiles.wx.begin(ff1), exv_distance_profiles.wx.end(ff1), sinqd_table->begin(q), 0.0);
                Iq[q-q0] -= 2*cx*this->free_params.cw*wx_sum*ff_ax_table.index(form_factor::water_bin, ff1).evaluate(q);
            }

            // water-water
            double ww_sum = std::inner_product(this->distance_profiles.ww.begin(), this->distance_profiles.ww.end(), sinqd_table->begin(q), 0.0);
            Iq[q-q0] += std::pow(this->free_params.cw, 2)*ww_sum*ff_aa_table.index(form_factor::water_bin, form_factor::water_bin).evaluate(q);
        }
        return ScatteringProfile(std::move(Iq), debye_axis);
    }
};

struct DebugCompositeDistanceHistogramFFGrid : public CompositeDistanceHistogramFFGrid {
    using CompositeDistanceHistogramFFGrid::CompositeDistanceHistogramFFGrid;
    using CompositeDistanceHistogramFFGrid::cache_get_intensity_profiles;
    DebugCompositeDistanceHistogramFFGrid(CompositeDistanceHistogramFFGrid&& other) : CompositeDistanceHistogramFFGrid(std::move(other)) {}

    inline static bool called_aa = false;
    inline static bool called_aw = false;
    inline static bool called_ww = false;
    inline static bool called_ax = false;
    inline static bool called_wx = false;
    inline static bool called_xx = false;
    inline static bool called_debye = false;
    ScatteringProfile get_profile_aa() const override  {
        called_aa = true;
        const auto& ff_table = get_ff_table();
        auto sinqd_table = sinc_table.get_sinc_table();
        Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
        unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

        std::vector<double> Iq(debye_axis.bins, 0);
        for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
            for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
                for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                    double aa_sum = std::inner_product(distance_profiles.aa.begin(ff1, ff2), distance_profiles.aa.end(ff1, ff2), sinqd_table->begin(q), 0.0);
                    Iq[q-q0] += aa_sum*ff_table.index(ff1, ff2).evaluate(q);
                }
            }
        }
        return ScatteringProfile(std::move(Iq), debye_axis);
    }

    ScatteringProfile get_profile_aw() const override  {
        called_aw = true;
        const auto& ff_table = get_ff_table();
        auto sinqd_table = sinc_table.get_sinc_table();
        Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
        unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

        std::vector<double> Iq(debye_axis.bins, 0);
        for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
            for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
                double aw_sum = std::inner_product(distance_profiles.aw.begin(ff1), distance_profiles.aw.end(ff1), sinqd_table->begin(q), 0.0);
                Iq[q-q0] += 2*free_params.cw*aw_sum*ff_table.index(ff1, form_factor::water_bin).evaluate(q);
            }
        }
        return ScatteringProfile(std::move(Iq), debye_axis);
    }

    ScatteringProfile get_profile_ww() const override  {
        called_ww = true;
        const auto& ff_table = get_ff_table();
        auto sinqd_table = sinc_table.get_sinc_table();
        Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
        unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

        std::vector<double> Iq(debye_axis.bins, 0);
        for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
            double ww_sum = std::inner_product(distance_profiles.ww.begin(), distance_profiles.ww.end(), sinqd_table->begin(q), 0.0);
            Iq[q-q0] += free_params.cw*free_params.cw*ww_sum*ff_table.index(form_factor::water_bin, form_factor::water_bin).evaluate(q);
        }
        return ScatteringProfile(std::move(Iq), debye_axis);
    }

    ScatteringProfile get_profile_ax() const override  {
        called_ax = true;
        const auto& ff_table = get_ff_table();
        auto sinqd_table = get_sinc_table_ax();
        Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
        unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

        std::vector<double> Iq(debye_axis.bins, 0);
        for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
            double cx = exv_factor(q);
            for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
                double ax_sum = std::inner_product(distance_profiles.aa.begin(ff1, form_factor::exv_bin), distance_profiles.aa.end(ff1, form_factor::exv_bin), sinqd_table->begin(q), 0.0);
                Iq[q-q0] += 2*cx*ax_sum*ff_table.index(ff1, form_factor::exv_bin).evaluate(q);
            }
        }
        return ScatteringProfile(std::move(Iq), debye_axis);
    }

    ScatteringProfile get_profile_wx() const override  {
        called_wx = true;
        const auto& ff_table = get_ff_table();
        auto sinqd_table = get_sinc_table_ax();
        Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
        unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

        std::vector<double> Iq(debye_axis.bins, 0);
        for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
            double cx = exv_factor(q);
            double ew_sum = std::inner_product(distance_profiles.aw.begin(form_factor::exv_bin), distance_profiles.aw.end(form_factor::exv_bin), sinqd_table->begin(q), 0.0);
            Iq[q-q0] += 2*cx*free_params.cw*ew_sum*ff_table.index(form_factor::exv_bin, form_factor::water_bin).evaluate(q);
        }
        return ScatteringProfile(std::move(Iq), debye_axis);
    }

    ScatteringProfile get_profile_xx() const override  {
        called_xx = true;
        const auto& ff_table = get_ff_table();
        auto sinqd_table = get_sinc_table_xx();
        Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
        unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

        std::vector<double> Iq(debye_axis.bins, 0);
        for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
            double cx = exv_factor(q);
            double xx_sum = std::inner_product(distance_profiles.aa.begin(form_factor::exv_bin, form_factor::exv_bin), distance_profiles.aa.end(form_factor::exv_bin, form_factor::exv_bin), sinqd_table->begin(q), 0.0);
            Iq[q-q0] += cx*cx*xx_sum*ff_table.index(form_factor::exv_bin, form_factor::exv_bin).evaluate(q);
        }
        return ScatteringProfile(std::move(Iq), debye_axis);
    }

    ScatteringProfile debye_transform() const override  {
        called_debye = true;
        const auto& ff_table = get_ff_table();
        auto sinqd_table_aa = sinc_table.get_sinc_table();
        auto sinqd_table_ax = get_sinc_table_ax();
        auto sinqd_table_xx = get_sinc_table_xx();

        // calculate the Debye scattering intensity
        Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
        unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

        std::vector<double> Iq(debye_axis.bins, 0);
        for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
            double cx = exv_factor(q);
            for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
                // atom-atom
                for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                    double aa_sum = std::inner_product(distance_profiles.aa.begin(ff1, ff2), distance_profiles.aa.end(ff1, ff2), sinqd_table_aa->begin(q), 0.0);
                    Iq[q-q0] += aa_sum*ff_table.index(ff1, ff2).evaluate(q);
                }

                // atom-exv
                double ax_sum = std::inner_product(distance_profiles.aa.begin(ff1, form_factor::exv_bin), distance_profiles.aa.end(ff1, form_factor::exv_bin), sinqd_table_ax->begin(q), 0.0);
                Iq[q-q0] -= 2*cx*ax_sum*ff_table.index(ff1, form_factor::exv_bin).evaluate(q);

                // atom-water
                double aw_sum = std::inner_product(distance_profiles.aw.begin(ff1), distance_profiles.aw.end(ff1), sinqd_table_aa->begin(q), 0.0);
                Iq[q-q0] += 2*free_params.cw*aw_sum*ff_table.index(ff1, form_factor::water_bin).evaluate(q);
            }

            // exv-exv
            double xx_sum = std::inner_product(distance_profiles.aa.begin(form_factor::exv_bin, form_factor::exv_bin), distance_profiles.aa.end(form_factor::exv_bin, form_factor::exv_bin), sinqd_table_xx->begin(q), 0.0);
            Iq[q-q0] += cx*cx*xx_sum*ff_table.index(form_factor::exv_bin, form_factor::exv_bin).evaluate(q);

            // exv-water
            double wx_sum = std::inner_product(distance_profiles.aw.begin(form_factor::exv_bin), distance_profiles.aw.end(form_factor::exv_bin), sinqd_table_ax->begin(q), 0.0);
            Iq[q-q0] -= 2*cx*free_params.cw*wx_sum*ff_table.index(form_factor::exv_bin, form_factor::water_bin).evaluate(q);

            // water-water
            double ww_sum = std::inner_product(distance_profiles.ww.begin(), distance_profiles.ww.end(), sinqd_table_aa->begin(q), 0.0);
            Iq[q-q0] += std::pow(free_params.cw, 2)*ww_sum*ff_table.index(form_factor::water_bin, form_factor::water_bin).evaluate(q);
        }
        return ScatteringProfile(std::move(Iq), debye_axis);
    }
};

struct DebugCompositeDistanceHistogramFFGridSurface : public CompositeDistanceHistogramFFGridSurface {
    using CompositeDistanceHistogramFFGridSurface::CompositeDistanceHistogramFFGridSurface;
    using CompositeDistanceHistogramFFGridSurface::cache_get_intensity_profiles;
    DebugCompositeDistanceHistogramFFGridSurface(CompositeDistanceHistogramFFGridSurface&& other) : CompositeDistanceHistogramFFGridSurface(std::move(other)) {}

    inline static bool called_aa = false;
    inline static bool called_aw = false;
    inline static bool called_ww = false;
    inline static bool called_ax = false;
    inline static bool called_wx = false;
    inline static bool called_xx = false;
    inline static bool called_debye = false;
    ScatteringProfile get_profile_aa() const override {
        called_aa = true;
        const auto& ff_table = get_ff_table();
        auto sinqd_table = sinc_table.get_sinc_table();
        Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
        unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

        std::vector<double> Iq(debye_axis.bins, 0);
        for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
            for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
                for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                    double aa_sum = std::inner_product(distance_profiles.aa.begin(ff1, ff2), distance_profiles.aa.end(ff1, ff2), sinqd_table->begin(q), 0.0);
                    Iq[q-q0] += aa_sum*ff_table.index(ff1, ff2).evaluate(q);
                }
            }
        }
        return ScatteringProfile(std::move(Iq), debye_axis);
    }

    ScatteringProfile get_profile_aw() const override {
        called_aw = true;
        const auto& ff_table = get_ff_table();
        auto sinqd_table = sinc_table.get_sinc_table();
        Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
        unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

        std::vector<double> Iq(debye_axis.bins, 0);
        for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
            for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
                double aw_sum = std::inner_product(distance_profiles.aw.begin(ff1), distance_profiles.aw.end(ff1), sinqd_table->begin(q), 0.0);
                Iq[q-q0] += 2*free_params.cw*aw_sum*ff_table.index(ff1, form_factor::water_bin).evaluate(q);
            }
        }
        return ScatteringProfile(std::move(Iq), debye_axis);
    }

    ScatteringProfile get_profile_ww() const override  {
        called_ww = true;
        const auto& ff_table = get_ff_table();
        auto sinqd_table = sinc_table.get_sinc_table();
        Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
        unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

        std::vector<double> Iq(debye_axis.bins, 0);
        for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
            double ww_sum = std::inner_product(distance_profiles.ww.begin(), distance_profiles.ww.end(), sinqd_table->begin(q), 0.0);
            Iq[q-q0] += free_params.cw*free_params.cw*ww_sum*ff_table.index(form_factor::water_bin, form_factor::water_bin).evaluate(q);
        }
        return ScatteringProfile(std::move(Iq), debye_axis);
    }

    ScatteringProfile get_profile_ax() const override  {
        called_ax = true;
        const auto& ff_table = get_ff_table();
        auto sinqd_table = get_sinc_table_ax();

        Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
        unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

        std::vector<double> Iq(debye_axis.bins, 0);
        for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
            auto ax = evaluate_ax_distance_profile(exv_factor(constants::axes::q_vals[q]));
            for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
                double ax_sum = std::inner_product(ax.begin(ff1), ax.end(ff1), sinqd_table->begin(q), 0.0);
                Iq[q-q0] += 2*ax_sum*ff_table.index(ff1, form_factor::exv_bin).evaluate(q);
            }
        }
        return ScatteringProfile(std::move(Iq), debye_axis);
    }

    ScatteringProfile get_profile_wx() const override  {
        called_wx = true;
        const auto& ff_table = get_ff_table();
        auto sinqd_table = get_sinc_table_ax();

        Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
        unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

        std::vector<double> Iq(debye_axis.bins, 0);
        for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
            auto wx = evaluate_wx_distance_profile(exv_factor(constants::axes::q_vals[q]));
            double ew_sum = std::inner_product(wx.begin(), wx.end(), sinqd_table->begin(q), 0.0);
            Iq[q-q0] += 2*free_params.cw*ew_sum*ff_table.index(form_factor::exv_bin, form_factor::water_bin).evaluate(q);
        }
        return ScatteringProfile(std::move(Iq), debye_axis);
    }

    ScatteringProfile get_profile_xx() const override  {
        called_xx = true;
        const auto& ff_table = get_ff_table();
        auto sinqd_table = get_sinc_table_xx();

        Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
        unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

        std::vector<double> Iq(debye_axis.bins, 0);
        for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
            auto xx = evaluate_xx_distance_profile(exv_factor(constants::axes::q_vals[q]));
            double xx_sum = std::inner_product(xx.begin(), xx.end(), sinqd_table->begin(q), 0.0);
            Iq[q-q0] += xx_sum*ff_table.index(form_factor::exv_bin, form_factor::exv_bin).evaluate(q);
        }
        return ScatteringProfile(std::move(Iq), debye_axis);
    }

    ScatteringProfile debye_transform() const override {
        called_debye = true;
        const auto& ff_table = get_ff_table();
        auto sinqd_table_aa = sinc_table.get_sinc_table();
        auto sinqd_table_ax = get_sinc_table_ax();
        auto sinqd_table_xx = get_sinc_table_xx();

        // calculate the Debye scattering intensity
        Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
        unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

        std::vector<double> Iq(debye_axis.bins, 0);
        for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
            double cx = exv_factor(constants::axes::q_vals[q]);
            auto xx = evaluate_xx_distance_profile(cx);
            auto wx = evaluate_wx_distance_profile(cx);
            auto ax = evaluate_ax_distance_profile(cx);

            for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
                // atom-atom
                for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                    double aa_sum = std::inner_product(distance_profiles.aa.begin(ff1, ff2), distance_profiles.aa.end(ff1, ff2), sinqd_table_aa->begin(q), 0.0);
                    Iq[q-q0] += aa_sum*ff_table.index(ff1, ff2).evaluate(q);
                }

                // atom-exv
                double ax_sum = std::inner_product(ax.begin(ff1), ax.end(ff1), sinqd_table_ax->begin(q), 0.0);
                Iq[q-q0] -= 2*ax_sum*ff_table.index(ff1, form_factor::exv_bin).evaluate(q);

                // atom-water
                double aw_sum = std::inner_product(distance_profiles.aw.begin(ff1), distance_profiles.aw.end(ff1), sinqd_table_aa->begin(q), 0.0);
                Iq[q-q0] += 2*free_params.cw*aw_sum*ff_table.index(ff1, form_factor::water_bin).evaluate(q);
            }

            // exv-exv
            double xx_sum = std::inner_product(xx.begin(), xx.end(), sinqd_table_xx->begin(q), 0.0);
            Iq[q-q0] += xx_sum*ff_table.index(form_factor::exv_bin, form_factor::exv_bin).evaluate(q);

            // exv-water
            double wx_sum = std::inner_product(wx.begin(), wx.end(), sinqd_table_ax->begin(q), 0.0);
            Iq[q-q0] -= 2*free_params.cw*wx_sum*ff_table.index(form_factor::exv_bin, form_factor::water_bin).evaluate(q);

            // water-water
            double ww_sum = std::inner_product(distance_profiles.ww.begin(), distance_profiles.ww.end(), sinqd_table_aa->begin(q), 0.0);
            Iq[q-q0] += std::pow(free_params.cw, 2)*ww_sum*ff_table.index(form_factor::water_bin, form_factor::water_bin).evaluate(q);
        }
        return ScatteringProfile(std::move(Iq), debye_axis);
    }
};

// Test the legacy single-threaded implementation of the Debye calculations against the new cached multi-threaded implementation
TEST_CASE("CompositeDistanceHistogramFFAvg: legacy comparison") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    std::string test_files = GENERATE("tests/files/2epe.pdb", "tests/files/c60.pdb", "tests/files/diamond.pdb");
    data::Molecule protein(test_files);

    SECTION(test_files + " full profile") {
        auto params = GENERATE(std::pair{1, 1}, std::pair{1, 2}, std::pair{2, 1}, std::pair{2, 2});

        auto hist_new = HistogramManagerMTFFAvg<false, false>(&protein).calculate_all();
        auto new_cast = static_cast<CompositeDistanceHistogramFFAvg*>(hist_new.get());
        new_cast->apply_water_scaling_factor(params.first);
        new_cast->apply_excluded_volume_scaling_factor(params.second);
        auto new_profile = hist_new->debye_transform();

        auto hist_old = DebugCompositeDistanceHistogramFFAvg(std::move(*new_cast));
        hist_old.apply_water_scaling_factor(params.first);
        hist_old.apply_excluded_volume_scaling_factor(params.second);
        auto old_profile = hist_old.debye_transform();

        SECTION("cw = " + std::to_string(params.first) + ", cx = " + std::to_string(params.second)) {
            CHECK(compare_hist(new_profile, old_profile));
            REQUIRE(DebugCompositeDistanceHistogramFFAvg::called_debye);
        }
    }

    SECTION(test_files + " individual profiles") {
        auto hist_new = HistogramManagerMTFFAvg<false, false>(&protein).calculate_all();
        auto new_cast = static_cast<CompositeDistanceHistogramFFAvg*>(hist_new.get());
        auto[_aa, _ax, _aw, _xx, _wx, _ww] = new_cast->cache_get_intensity_profiles();

        auto hist_old = DebugCompositeDistanceHistogramFFAvg(std::move(*new_cast));
        auto ww = hist_old.get_profile_ww();
        auto wx = hist_old.get_profile_wx();
        auto xx = hist_old.get_profile_xx();
        auto aw = hist_old.get_profile_aw();
        auto ax = hist_old.get_profile_ax();
        auto aa = hist_old.get_profile_aa();

        CHECK(compare_hist(ww, _ww));
        CHECK(compare_hist(wx, _wx));
        CHECK(compare_hist(xx, _xx));
        CHECK(compare_hist(aw, _aw));
        CHECK(compare_hist(ax, _ax));
        CHECK(compare_hist(aa, _aa));

        REQUIRE(DebugCompositeDistanceHistogramFFAvg::called_aa);
        REQUIRE(DebugCompositeDistanceHistogramFFAvg::called_aw);
        REQUIRE(DebugCompositeDistanceHistogramFFAvg::called_ww);
        REQUIRE(DebugCompositeDistanceHistogramFFAvg::called_ax);
        REQUIRE(DebugCompositeDistanceHistogramFFAvg::called_wx);
        REQUIRE(DebugCompositeDistanceHistogramFFAvg::called_xx);
    }
}

// Test the legacy single-threaded implementation of the Debye calculations against the new cached multi-threaded implementation
TEST_CASE("CompositeDistanceHistogramFFExplicit: legacy comparison") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    std::string test_files = GENERATE("tests/files/2epe.pdb", "tests/files/c60.pdb", "tests/files/diamond.pdb");
    data::Molecule protein(test_files);

    SECTION(test_files + " full profile") {
        auto params = GENERATE(std::pair{1, 1}, std::pair{1, 2}, std::pair{2, 1}, std::pair{2, 2});

        auto hist_new = HistogramManagerMTFFExplicit<false, false>(&protein).calculate_all();
        auto new_cast = static_cast<CompositeDistanceHistogramFFExplicit*>(hist_new.get());
        new_cast->apply_water_scaling_factor(params.first);
        new_cast->apply_excluded_volume_scaling_factor(params.second);
        auto new_profile = hist_new->debye_transform();

        auto hist_old = DebugCompositeDistanceHistogramFFExplicit(std::move(*new_cast));
        hist_old.apply_water_scaling_factor(params.first);
        hist_old.apply_excluded_volume_scaling_factor(params.second);
        auto old_profile = hist_old.debye_transform();

        SECTION("cw = " + std::to_string(params.first) + ", cx = " + std::to_string(params.second)) {
            CHECK(compare_hist(new_profile, old_profile));
            REQUIRE(DebugCompositeDistanceHistogramFFExplicit::called_debye);
        }
    }

    SECTION(test_files + " individual profiles") {
        auto hist_new = HistogramManagerMTFFExplicit<false, false>(&protein).calculate_all();
        auto new_cast = static_cast<CompositeDistanceHistogramFFExplicit*>(hist_new.get());
        auto[_aa, _ax, _aw, _xx, _wx, _ww] = new_cast->cache_get_intensity_profiles();

        auto hist_old = DebugCompositeDistanceHistogramFFExplicit(std::move(*new_cast));
        auto ww = hist_old.get_profile_ww();
        auto wx = hist_old.get_profile_wx();
        auto xx = hist_old.get_profile_xx();
        auto aw = hist_old.get_profile_aw();
        auto ax = hist_old.get_profile_ax();
        auto aa = hist_old.get_profile_aa();

        CHECK(compare_hist(ww, _ww));
        CHECK(compare_hist(wx, _wx));
        CHECK(compare_hist(xx, _xx));
        CHECK(compare_hist(aw, _aw));
        CHECK(compare_hist(ax, _ax));
        CHECK(compare_hist(aa, _aa));

        REQUIRE(DebugCompositeDistanceHistogramFFExplicit::called_aa);
        REQUIRE(DebugCompositeDistanceHistogramFFExplicit::called_aw);
        REQUIRE(DebugCompositeDistanceHistogramFFExplicit::called_ww);
        REQUIRE(DebugCompositeDistanceHistogramFFExplicit::called_ax);
        REQUIRE(DebugCompositeDistanceHistogramFFExplicit::called_wx);
        REQUIRE(DebugCompositeDistanceHistogramFFExplicit::called_xx);
    }
}

// Test the legacy single-threaded implementation of the Debye calculations against the new cached multi-threaded implementation
TEST_CASE("CompositeDistanceHistogramFFGrid: legacy comparison") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    std::string test_files = GENERATE("tests/files/2epe.pdb", "tests/files/c60.pdb", "tests/files/diamond.pdb");
    data::Molecule protein(test_files);

    SECTION(test_files + " full profile") {
        auto params = GENERATE(std::pair{1, 1}, std::pair{1, 2}, std::pair{2, 1}, std::pair{2, 2});

        auto hist_new = HistogramManagerMTFFGrid<false>(&protein).calculate_all();
        auto new_cast = static_cast<CompositeDistanceHistogramFFGrid*>(hist_new.get());
        new_cast->apply_water_scaling_factor(params.first);
        new_cast->apply_excluded_volume_scaling_factor(params.second);
        auto new_profile = hist_new->debye_transform();

        auto hist_old = DebugCompositeDistanceHistogramFFGrid(std::move(*new_cast));
        hist_old.apply_water_scaling_factor(params.first);
        hist_old.apply_excluded_volume_scaling_factor(params.second);
        auto old_profile = hist_old.debye_transform();

        SECTION("cw = " + std::to_string(params.first) + ", cx = " + std::to_string(params.second)) {
            CHECK(compare_hist(new_profile, old_profile));
            REQUIRE(DebugCompositeDistanceHistogramFFGrid::called_debye);
        }
    }

    SECTION(test_files + " individual profiles") {
        auto hist_new = HistogramManagerMTFFGrid<false>(&protein).calculate_all();
        auto new_cast = static_cast<CompositeDistanceHistogramFFGrid*>(hist_new.get());
        auto[_aa, _ax, _aw, _xx, _wx, _ww] = new_cast->cache_get_intensity_profiles();

        auto hist_old = DebugCompositeDistanceHistogramFFGrid(std::move(*new_cast));
        auto ww = hist_old.get_profile_ww();
        auto wx = hist_old.get_profile_wx();
        auto xx = hist_old.get_profile_xx();
        auto aw = hist_old.get_profile_aw();
        auto ax = hist_old.get_profile_ax();
        auto aa = hist_old.get_profile_aa();

        CHECK(compare_hist(ww, _ww));
        CHECK(compare_hist(wx, _wx));
        CHECK(compare_hist(xx, _xx));
        CHECK(compare_hist(aw, _aw));
        CHECK(compare_hist(ax, _ax));
        CHECK(compare_hist(aa, _aa));

        REQUIRE(DebugCompositeDistanceHistogramFFGrid::called_aa);
        REQUIRE(DebugCompositeDistanceHistogramFFGrid::called_aw);
        REQUIRE(DebugCompositeDistanceHistogramFFGrid::called_ww);
        REQUIRE(DebugCompositeDistanceHistogramFFGrid::called_ax);
        REQUIRE(DebugCompositeDistanceHistogramFFGrid::called_wx);
        REQUIRE(DebugCompositeDistanceHistogramFFGrid::called_xx);
    }
}

// Test the legacy single-threaded implementation of the Debye calculations against the new cached multi-threaded implementation
TEST_CASE("CompositeDistanceHistogramFFGridSurface: legacy comparison") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    std::string test_files = GENERATE("tests/files/2epe.pdb", "tests/files/c60.pdb", "tests/files/diamond.pdb");
    data::Molecule protein(test_files);

    SECTION(test_files + " full profile") {
        auto params = GENERATE(std::pair{1, 1}, std::pair{1, 2}, std::pair{2, 1}, std::pair{2, 2});

        auto hist_new = HistogramManagerMTFFGridSurface<false>(&protein).calculate_all();
        auto new_cast = static_cast<CompositeDistanceHistogramFFGridSurface*>(hist_new.get());
        new_cast->apply_water_scaling_factor(params.first);
        new_cast->apply_excluded_volume_scaling_factor(params.second);
        auto new_profile = hist_new->debye_transform();

        auto hist_old = DebugCompositeDistanceHistogramFFGridSurface(std::move(*new_cast));
        hist_old.apply_water_scaling_factor(params.first);
        hist_old.apply_excluded_volume_scaling_factor(params.second);
        auto old_profile = hist_old.debye_transform();

        SECTION("cw = " + std::to_string(params.first) + ", cx = " + std::to_string(params.second)) {
            CHECK(compare_hist(new_profile, old_profile));
            REQUIRE(DebugCompositeDistanceHistogramFFGridSurface::called_debye);
        }
    }

    SECTION(test_files + " individual profiles") {
        auto hist_new = HistogramManagerMTFFGridSurface<false>(&protein).calculate_all();
        auto new_cast = static_cast<CompositeDistanceHistogramFFGridSurface*>(hist_new.get());
        auto[_aa, _ax, _aw, _xx, _wx, _ww] = new_cast->cache_get_intensity_profiles();

        auto hist_old = DebugCompositeDistanceHistogramFFGridSurface(std::move(*new_cast));
        auto ww = hist_old.get_profile_ww();
        auto wx = hist_old.get_profile_wx();
        auto xx = hist_old.get_profile_xx();
        auto aw = hist_old.get_profile_aw();
        auto ax = hist_old.get_profile_ax();
        auto aa = hist_old.get_profile_aa();

        CHECK(compare_hist(ww, _ww));
        CHECK(compare_hist(wx, _wx));
        CHECK(compare_hist(xx, _xx));
        CHECK(compare_hist(aw, _aw));
        CHECK(compare_hist(ax, _ax));
        CHECK(compare_hist(aa, _aa));

        REQUIRE(DebugCompositeDistanceHistogramFFGridSurface::called_aa);
        REQUIRE(DebugCompositeDistanceHistogramFFGridSurface::called_aw);
        REQUIRE(DebugCompositeDistanceHistogramFFGridSurface::called_ww);
        REQUIRE(DebugCompositeDistanceHistogramFFGridSurface::called_ax);
        REQUIRE(DebugCompositeDistanceHistogramFFGridSurface::called_wx);
        REQUIRE(DebugCompositeDistanceHistogramFFGridSurface::called_xx);
    }
}