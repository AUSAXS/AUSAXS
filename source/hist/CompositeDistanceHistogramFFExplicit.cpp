#include <hist/CompositeDistanceHistogramFFExplicit.h>
#include <hist/Histogram.h>
#include <table/ArrayDebyeTable.h>
#include <form_factor/FormFactor.h>
#include <form_factor/PrecalculatedFormFactorProduct.h>
#include <form_factor/PrecalculatedExvFormFactorProduct.h>
#include <settings/HistogramSettings.h>

using namespace hist;

CompositeDistanceHistogramFFExplicit::CompositeDistanceHistogramFFExplicit() = default;

CompositeDistanceHistogramFFExplicit::CompositeDistanceHistogramFFExplicit(container::Container3D<double>&& p_aa, container::Container3D<double>&& p_ax, container::Container3D<double>&& p_xx,
                                                                           container::Container2D<double>&& p_wa, container::Container2D<double>&& p_wx, container::Container1D<double>&& p_ww,
                                                                           std::vector<double>&& p_tot, const Axis& axis)
    : CompositeDistanceHistogramFF(std::move(p_aa), std::move(p_wa), std::move(p_ww), std::move(p_tot), axis), cp_ax(std::move(p_ax)), cp_xx(std::move(p_xx)), cp_wx(std::move(p_wx)) {}

CompositeDistanceHistogramFFExplicit::~CompositeDistanceHistogramFFExplicit() = default;

ScatteringProfile CompositeDistanceHistogramFFExplicit::debye_transform() const {
    const auto& ff_aa_table = form_factor::storage::get_precalculated_form_factor_table();
    const auto& ff_ax_table = form_factor::storage::cross::get_precalculated_form_factor_table();
    const auto& ff_xx_table = form_factor::storage::exv::get_precalculated_form_factor_table();
    const auto& sinqd_table = table::ArrayDebyeTable::get_default_table();

    // calculate the Debye scattering intensity
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    std::vector<double> q_axis = debye_axis.as_vector();

    unsigned int ff_w_index = static_cast<int>(form_factor::form_factor_t::OH);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = ff1; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                // atom-atom
                double aa_sum = std::inner_product(cp_aa.begin(ff1, ff2), cp_aa.end(ff1, ff2), sinqd_table.begin(q), 0.0);
                Iq[q] += aa_sum*ff_aa_table.index(ff1, ff2).evaluate(q);

                // atom-exv
                double ax_sum = std::inner_product(cp_ax.begin(ff1, ff2), cp_ax.end(ff1, ff2), sinqd_table.begin(q), 0.0);
                Iq[q] += ax_sum*ff_ax_table.index(ff1, ff2).evaluate(q);

                // exv-exv
                double exv_exv_sum = std::inner_product(cp_xx.begin(ff1, ff2), cp_xx.end(ff1, ff2), sinqd_table.begin(q), 0.0);
                Iq[q] += exv_exv_sum*ff_xx_table.index(ff1, ff2).evaluate(q);
            }

            // atom-water
            double aw_sum = std::inner_product(cp_wa.begin(ff1), cp_wa.end(ff1), sinqd_table.begin(q), 0.0);
            Iq[q] += 2*w_scaling*aw_sum*ff_aa_table.index(ff1, ff_w_index).evaluate(q);

            // exv-water
            double wx_sum = std::inner_product(cp_wx.begin(ff1), cp_wx.end(ff1), sinqd_table.begin(q), 0.0);
            Iq[q] += w_scaling*wx_sum*ff_ax_table.index(ff1, ff_w_index).evaluate(q);
        }

        // water-water
        double ww_sum = std::inner_product(cp_ww.begin(), cp_ww.end(), sinqd_table.begin(q), 0.0);
        Iq[q] += w_scaling*w_scaling*ww_sum*ff_aa_table.index(ff_w_index, ff_w_index).evaluate(q);
    }
    return ScatteringProfile(Iq, debye_axis);
}

const std::vector<double>& CompositeDistanceHistogramFF::get_counts() const {
    p = std::vector<double>(axis.bins, 0);
    auto& p_pp = get_pp_counts();
    auto& p_hp = get_hp_counts();
    auto& p_hh = get_hh_counts();
    for (unsigned int i = 0; i < axis.bins; ++i) {
        p[i] = p_pp[i] + 2*w_scaling*p_hp[i] + w_scaling*w_scaling*p_hh[i];
    }
    return p.data;
}

const std::vector<double>& CompositeDistanceHistogramFF::get_pp_counts() const {
    p_pp = std::vector<double>(axis.bins, 0);
    for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
        for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
            std::transform(p_pp.begin(), p_pp.end(), cp_aa.begin(ff1, ff2), p_pp.begin(), std::plus<double>());
        }
    }
    return p_pp;
}

const std::vector<double>& CompositeDistanceHistogramFF::get_hp_counts() const {
    p_hp = std::vector<double>(axis.bins, 0);
    for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
        std::transform(p_hp.begin(), p_hp.end(), cp_hp.begin(ff1), p_hp.begin(), std::plus<double>());
    }
    return p_hp;
}

const std::vector<double>& CompositeDistanceHistogramFF::get_hh_counts() const {
    p_hh = std::vector<double>(axis.bins, 0);
    std::transform(p_hh.begin(), p_hh.end(), cp_hh.begin(), p_hh.begin(), std::plus<double>());
    return p_hh;
}

void CompositeDistanceHistogramFF::apply_water_scaling_factor(double k) {
    w_scaling = k;
}

void CompositeDistanceHistogramFF::apply_excluded_volume_scaling_factor(double k) {
    exv_scaling = k;
}