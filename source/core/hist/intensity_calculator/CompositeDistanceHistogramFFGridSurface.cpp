/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hist/intensity_calculator/CompositeDistanceHistogramFFGridSurface.h>
#include <form_factor/PrecalculatedFormFactorProduct.h>
#include <form_factor/ExvFormFactor.h>
#include <table/ArrayDebyeTable.h>
#include <settings/GridSettings.h>
#include <settings/HistogramSettings.h>
#include <utility/Exceptions.h>
#include <dataset/SimpleDataset.h>

using namespace hist;
using namespace form_factor;

CompositeDistanceHistogramFFGridSurface::CompositeDistanceHistogramFFGridSurface(CompositeDistanceHistogramFFGridSurface&&) noexcept = default;
CompositeDistanceHistogramFFGridSurface& CompositeDistanceHistogramFFGridSurface::operator=(CompositeDistanceHistogramFFGridSurface&&) noexcept = default;
CompositeDistanceHistogramFFGridSurface::~CompositeDistanceHistogramFFGridSurface() = default;

CompositeDistanceHistogramFFGridSurface::XXContainer CompositeDistanceHistogramFFGridSurface::XXContainer::operator+=(const CompositeDistanceHistogramFFGridSurface::XXContainer& other) {
    std::transform(interior.begin(), interior.end(), other.interior.begin(), interior.begin(), std::plus<>());
    std::transform(surface.begin(), surface.end(), other.surface.begin(), surface.begin(), std::plus<>());
    std::transform(cross.begin(), cross.end(), other.cross.begin(), cross.begin(), std::plus<>());
    return *this;
}

CompositeDistanceHistogramFFGridSurface::AXContainer CompositeDistanceHistogramFFGridSurface::AXContainer::operator+=(const CompositeDistanceHistogramFFGridSurface::AXContainer& other) {
    std::transform(interior.begin(), interior.end(), other.interior.begin(), interior.begin(), [](auto& a, const auto& b) {return a+b;});
    std::transform(surface.begin(), surface.end(), other.surface.begin(), surface.begin(), [](auto& a, const auto& b) {return a+b;});
    return *this;
}

CompositeDistanceHistogramFFGridSurface::WXContainer CompositeDistanceHistogramFFGridSurface::WXContainer::operator+=(const CompositeDistanceHistogramFFGridSurface::WXContainer& other) {
    std::transform(interior.begin(), interior.end(), other.interior.begin(), interior.begin(), std::plus<>());
    std::transform(surface.begin(), surface.end(), other.surface.begin(), surface.begin(), std::plus<>());
    return *this;
}

CompositeDistanceHistogramFFGridSurface::CompositeDistanceHistogramFFGridSurface(
    hist::Distribution3D&& p_aa, 
    hist::Distribution2D&& p_aw, 
    hist::Distribution1D&& p_ww, 
    XXContainer&& xx,
    AXContainer&& ax,
    WXContainer&& wx,
    hist::WeightedDistribution1D&& p_tot_aa,
    hist::WeightedDistribution1D&& p_tot_ax,
    hist::WeightedDistribution1D&& p_tot_xx
) : hist::CompositeDistanceHistogramFFAvg(std::move(p_aa), std::move(p_aw), std::move(p_ww), std::move(p_tot_aa)),
    exv_distance_profiles{hist::Distribution1D(std::move(xx.interior)), hist::Distribution1D(std::move(xx.surface)), hist::Distribution1D(std::move(xx.cross)), 
                          hist::Distribution1D(std::move(wx.interior)), hist::Distribution1D(std::move(wx.surface)), hist::Distribution2D(std::move(ax.interior)), 
                          hist::Distribution2D(std::move(ax.surface))}
{
    initialize(p_tot_ax.get_weighted_axis(), p_tot_xx.get_weighted_axis());
}

void CompositeDistanceHistogramFFGridSurface::regenerate_table() {ff_table = CompositeDistanceHistogramFFGrid::generate_table();}

form_factor::storage::atomic::table_t CompositeDistanceHistogramFFGridSurface::ff_table = CompositeDistanceHistogramFFGrid::generate_table();

Limit CompositeDistanceHistogramFFGridSurface::get_excluded_volume_scaling_factor_limits() const {
    return {0, 2};
}

double CompositeDistanceHistogramFFGridSurface::exv_factor(double q, double cx) {
    constexpr double rm2 = constants::radius::average_atomic_radius*constants::radius::average_atomic_radius;
    return std::pow(cx, 3)*std::exp(-rm2*(std::pow(cx, 2) - 1)*q*q/4);
}

double CompositeDistanceHistogramFFGridSurface::exv_factor(double q) const {
    return exv_factor(q, free_params.cx);
}

hist::Distribution1D CompositeDistanceHistogramFFGridSurface::evaluate_xx_profile(double cx) const {
    hist::Distribution1D xx = exv_distance_profiles.xx_i;
    std::transform(xx.begin(), xx.end(), exv_distance_profiles.xx_s.begin(), xx.begin(), [cx] (double a, double b) {return a + std::pow(cx, 2)*b;});
    std::transform(xx.begin(), xx.end(), exv_distance_profiles.xx_c.begin(), xx.begin(), [cx] (double a, double b) {return a + cx*b;});
    return xx;
}

hist::Distribution1D CompositeDistanceHistogramFFGridSurface::evaluate_wx_profile(double cx) const {
    hist::Distribution1D wx = exv_distance_profiles.wx_i;
    std::transform(wx.begin(), wx.end(), exv_distance_profiles.wx_s.begin(), wx.begin(), [cx] (double a, double b) {return a + cx*b;});
    return wx;
}

hist::Distribution2D CompositeDistanceHistogramFFGridSurface::evaluate_ax_profile(double cx) const {
    hist::Distribution2D ax = exv_distance_profiles.ax_i;
    std::transform(ax.begin(), ax.end(), exv_distance_profiles.ax_s.begin(), ax.begin(), [cx] (double a, double b) {return a + cx*b;});
    return ax;
}

ScatteringProfile CompositeDistanceHistogramFFGridSurface::debye_transform() const {
    const auto& ff_table = get_ff_table();
    auto sinqd_table_aa = get_sinc_table();
    auto sinqd_table_ax = get_sinc_table_ax();
    auto sinqd_table_xx = get_sinc_table_xx();

    // calculate the Debye scattering intensity
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        double cx = exv_factor(constants::axes::q_vals[q]);
        auto xx = evaluate_xx_profile(cx);
        auto wx = evaluate_wx_profile(cx);
        auto ax = evaluate_ax_profile(cx);

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

SimpleDataset CompositeDistanceHistogramFFGridSurface::debye_transform(const std::vector<double>&) const {
    throw except::not_implemented("CompositeDistanceHistogramFFGrid::debye_transform(const std::vector<double>& q) const");
}

ScatteringProfile CompositeDistanceHistogramFFGridSurface::get_profile_ax() const {
    const auto& ff_table = get_ff_table();
    auto sinqd_table = get_sinc_table_ax();

    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        auto ax = evaluate_ax_profile(exv_factor(constants::axes::q_vals[q]));
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            double ax_sum = std::inner_product(ax.begin(ff1), ax.end(ff1), sinqd_table->begin(q), 0.0);
            Iq[q-q0] += 2*ax_sum*ff_table.index(ff1, form_factor::exv_bin).evaluate(q);
        }
    }
    return ScatteringProfile(std::move(Iq), debye_axis);
}

ScatteringProfile CompositeDistanceHistogramFFGridSurface::get_profile_wx() const {
    const auto& ff_table = get_ff_table();
    auto sinqd_table = get_sinc_table_ax();

    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        auto wx = evaluate_wx_profile(exv_factor(constants::axes::q_vals[q]));
        double ew_sum = std::inner_product(wx.begin(), wx.end(), sinqd_table->begin(q), 0.0);
        Iq[q-q0] += 2*free_params.cw*ew_sum*ff_table.index(form_factor::exv_bin, form_factor::water_bin).evaluate(q);
    }
    return ScatteringProfile(std::move(Iq), debye_axis);
}

ScatteringProfile CompositeDistanceHistogramFFGridSurface::get_profile_xx() const {
    const auto& ff_table = get_ff_table();
    auto sinqd_table = get_sinc_table_xx();

    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        auto xx = evaluate_xx_profile(exv_factor(constants::axes::q_vals[q]));
        double xx_sum = std::inner_product(xx.begin(), xx.end(), sinqd_table->begin(q), 0.0);
        Iq[q-q0] += xx_sum*ff_table.index(form_factor::exv_bin, form_factor::exv_bin).evaluate(q);
    }
    return ScatteringProfile(std::move(Iq), debye_axis);
}

observer_ptr<const table::DebyeTable> CompositeDistanceHistogramFFGridSurface::get_sinc_table_ax() const {
    if (use_weighted_table) {return sinc_tables.ax.get();}
    return &table::ArrayDebyeTable::get_default_table();
}

observer_ptr<const table::DebyeTable> CompositeDistanceHistogramFFGridSurface::get_sinc_table_xx() const {
    if (use_weighted_table) {return sinc_tables.xx.get();}
    return &table::ArrayDebyeTable::get_default_table();
}

void CompositeDistanceHistogramFFGridSurface::initialize(std::vector<double>&& d_axis_ax, std::vector<double>&& d_axis_xx) {
    this->distance_axes = {.xx=std::move(d_axis_xx), .ax=std::move(d_axis_ax)};
    sinc_tables = {.xx=std::make_unique<table::VectorDebyeTable>(this->distance_axes.xx), .ax=std::make_unique<table::VectorDebyeTable>(this->distance_axes.ax)};

    // fix the aa counts to also contain the exv contributions
    auto xx = evaluate_xx_profile(1);
    auto wx = evaluate_wx_profile(1);
    auto ax = evaluate_ax_profile(1);

    auto& aa = CompositeDistanceHistogramFFAvgBase::get_aa_counts_ff();
    for (unsigned int ff = 0; ff < form_factor::get_count_without_excluded_volume(); ++ff) {
        std::transform(aa.begin(ff, form_factor::exv_bin), aa.end(ff, form_factor::exv_bin), ax.begin(ff), aa.begin(ff, form_factor::exv_bin), std::plus<double>());
    }
    std::transform(aa.begin(form_factor::exv_bin, form_factor::exv_bin), aa.end(form_factor::exv_bin, form_factor::exv_bin), xx.begin(), aa.begin(form_factor::exv_bin, form_factor::exv_bin), std::plus<double>());

    auto& aw = CompositeDistanceHistogramFFAvgBase::get_aw_counts_ff();
    std::transform(aw.begin(form_factor::exv_bin), aw.end(form_factor::exv_bin), wx.begin(), aw.begin(form_factor::exv_bin), std::plus<double>());
}