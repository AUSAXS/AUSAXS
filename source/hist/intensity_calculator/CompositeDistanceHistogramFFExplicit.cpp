/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hist/intensity_calculator/CompositeDistanceHistogramFFExplicit.h>
#include <hist/Histogram.h>
#include <table/ArrayDebyeTable.h>
#include <form_factor/FormFactor.h>
#include <form_factor/PrecalculatedFormFactorProduct.h>
#include <form_factor/PrecalculatedExvFormFactorProduct.h>
#include <settings/HistogramSettings.h>
#include <constants/ConstantsMath.h>
#include <math/ConstexprMath.h>

using namespace hist;

CompositeDistanceHistogramFFExplicit::CompositeDistanceHistogramFFExplicit() = default;

CompositeDistanceHistogramFFExplicit::CompositeDistanceHistogramFFExplicit(
    hist::Distribution3D&& p_aa, 
    hist::Distribution3D&& p_ax, 
    hist::Distribution3D&& p_xx, 
    hist::Distribution2D&& p_aw, 
    hist::Distribution2D&& p_wx, 
    hist::Distribution1D&& p_ww,
    hist::Distribution1D&& p_tot
) : CompositeDistanceHistogramFFAvg(std::move(p_aa), std::move(p_aw), std::move(p_ww), std::move(p_tot)), cp_ax(std::move(p_ax)), cp_xx(std::move(p_xx)), cp_wx(std::move(p_wx)) {}

CompositeDistanceHistogramFFExplicit::CompositeDistanceHistogramFFExplicit(
    hist::Distribution3D&& p_aa, 
    hist::Distribution3D&& p_ax, 
    hist::Distribution3D&& p_xx, 
    hist::Distribution2D&& p_aw, 
    hist::Distribution2D&& p_wx, 
    hist::Distribution1D&& p_ww, 
    hist::WeightedDistribution1D&& p_tot
) : CompositeDistanceHistogramFFAvg(std::move(p_aa), std::move(p_aw), std::move(p_ww), std::move(p_tot)), cp_ax(std::move(p_ax)), cp_xx(std::move(p_xx)), cp_wx(std::move(p_wx)) {}

CompositeDistanceHistogramFFExplicit::~CompositeDistanceHistogramFFExplicit() = default;

double CompositeDistanceHistogramFFExplicit::exv_factor(double q) const {
    // G(q) factor from CRYSOL
    constexpr double rm = 1.62;
    constexpr double c = constexpr_math::pow(4*constants::pi/3, 3./2)*constants::pi*rm*rm*constants::form_factor::s_to_q_factor;
    return std::pow(cx, 3)*std::exp(-c*(cx*cx - 1)*q*q);
}

// static unsigned int qcheck = 26;
// ScatteringProfile CompositeDistanceHistogramFFExplicit::debye_transform() const {
//     const auto& ff_aa_table = form_factor::storage::get_precalculated_form_factor_table();
//     const auto& ff_ax_table = form_factor::storage::cross::get_precalculated_form_factor_table();
//     const auto& ff_xx_table = form_factor::storage::exv::get_precalculated_form_factor_table();
//     const auto& sinqd_table = get_sinc_table();

//     // calculate the Debye scattering intensity
//     Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
//     unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

//     std::vector<double> Iq(debye_axis.bins, 0);
//     unsigned int ff_w_index = static_cast<int>(form_factor::form_factor_t::OH);
//     for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
//         double Gq = G_factor(constants::axes::q_vals[q]);
//         for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
//             for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
//                 // atom-atom
//                 double aa_sum = std::inner_product(cp_aa.begin(ff1, ff2), cp_aa.end(ff1, ff2), sinqd_table->begin(q), 0.0);
//                 Iq[q-q0] += aa_sum*ff_aa_table.index(ff1, ff2).evaluate(q);
//             }
//         }
//     }
//     std::cout << "(aa) I(" + std::to_string(qcheck) + ") = " << Iq[qcheck] << std::endl;
//     double tmp = Iq[qcheck];

//     for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
//         double Gq = G_factor(constants::axes::q_vals[q]);
//         for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
//             for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
//                 // atom-exv
//                 double ax_sum = std::inner_product(cp_ax.begin(ff1, ff2), cp_ax.end(ff1, ff2), sinqd_table->begin(q), 0.0);
//                 Iq[q-q0] -= 2*Gq*ax_sum*ff_ax_table.index(ff1, ff2).evaluate(q);
//             }
//         }
//     }
//     std::cout << "(ax) I(" + std::to_string(qcheck) + ") = " << Iq[qcheck] - tmp << std::endl;
//     tmp = Iq[qcheck];

//     for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
//         double Gq = G_factor(constants::axes::q_vals[q]);
//         for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
//             for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
//                 // exv-exv
//                 double xx_sum = std::inner_product(cp_xx.begin(ff1, ff2), cp_xx.end(ff1, ff2), sinqd_table->begin(q), 0.0);
//                 Iq[q-q0] += Gq*Gq*xx_sum*ff_xx_table.index(ff1, ff2).evaluate(q);
//             }
//         }
//     }
//     std::cout << "(xx) I(" + std::to_string(qcheck) + ") = " << Iq[qcheck] - tmp << std::endl;
//     tmp = Iq[qcheck];
    
//     for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
//         double Gq = G_factor(constants::axes::q_vals[q]);
//         for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
//             // atom-water
//             double aw_sum = std::inner_product(cp_aw.begin(ff1), cp_aw.end(ff1), sinqd_table->begin(q), 0.0);
//             Iq[q-q0] += 2*cw*aw_sum*ff_aa_table.index(ff1, ff_w_index).evaluate(q);
//         }
//     }
//     std::cout << "(aw) I(" + std::to_string(qcheck) + ") = " << Iq[qcheck] - tmp << std::endl;
//     tmp = Iq[qcheck];

//     for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
//         double Gq = G_factor(constants::axes::q_vals[q]);
//         for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
//             // exv-water
//             double wx_sum = std::inner_product(cp_wx.begin(ff1), cp_wx.end(ff1), sinqd_table->begin(q), 0.0);
//             Iq[q-q0] -= 2*Gq*cw*wx_sum*ff_ax_table.index(ff_w_index, ff1).evaluate(q);
//         }
//     }
//     std::cout << "(wx) I(" + std::to_string(qcheck) + ") = " << Iq[qcheck] - tmp << std::endl;
//     tmp = Iq[qcheck];

//     for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
//         double Gq = G_factor(constants::axes::q_vals[q]);
//         // water-water
//         double ww_sum = std::inner_product(cp_ww.begin(), cp_ww.end(), sinqd_table->begin(q), 0.0);
//         Iq[q-q0] += cw*cw*ww_sum*ff_aa_table.index(ff_w_index, ff_w_index).evaluate(q);
//     }
//     std::cout << "(ww) I(" + std::to_string(qcheck) + ") = " << Iq[qcheck] - tmp << std::endl;
//     tmp = Iq[qcheck];
//     return ScatteringProfile(Iq, debye_axis);
// }

ScatteringProfile CompositeDistanceHistogramFFExplicit::debye_transform() const {
    const auto& ff_aa_table = form_factor::storage::atomic::get_precalculated_form_factor_table();
    const auto& ff_ax_table = form_factor::storage::cross::get_precalculated_form_factor_table();
    const auto& ff_xx_table = form_factor::storage::exv::get_precalculated_form_factor_table();
    auto sinqd_table = get_sinc_table();

    // calculate the Debye scattering intensity
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    unsigned int ff_w_index = static_cast<int>(form_factor::form_factor_t::OH);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        double cx = exv_factor(constants::axes::q_vals[q]);
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                // atom-atom
                double aa_sum = std::inner_product(cp_aa.begin(ff1, ff2), cp_aa.end(ff1, ff2), sinqd_table->begin(q), 0.0);
                Iq[q-q0] += aa_sum*ff_aa_table.index(ff1, ff2).evaluate(q);

                // atom-exv
                double ax_sum = std::inner_product(cp_ax.begin(ff1, ff2), cp_ax.end(ff1, ff2), sinqd_table->begin(q), 0.0);
                Iq[q-q0] -= 2*cx*ax_sum*ff_ax_table.index(ff1, ff2).evaluate(q);

                // exv-exv
                double xx_sum = std::inner_product(cp_xx.begin(ff1, ff2), cp_xx.end(ff1, ff2), sinqd_table->begin(q), 0.0);
                Iq[q-q0] += cx*cx*xx_sum*ff_xx_table.index(ff1, ff2).evaluate(q);
            }

            // atom-water
            double aw_sum = std::inner_product(cp_aw.begin(ff1), cp_aw.end(ff1), sinqd_table->begin(q), 0.0);
            Iq[q-q0] += 2*cw*aw_sum*ff_aa_table.index(ff1, ff_w_index).evaluate(q);

            // exv-water
            double wx_sum = std::inner_product(cp_wx.begin(ff1), cp_wx.end(ff1), sinqd_table->begin(q), 0.0);
            Iq[q-q0] -= 2*cx*cw*wx_sum*ff_ax_table.index(ff_w_index, ff1).evaluate(q);
        }

        // water-water
        double ww_sum = std::inner_product(cp_ww.begin(), cp_ww.end(), sinqd_table->begin(q), 0.0);
        Iq[q-q0] += cw*cw*ww_sum*ff_aa_table.index(ff_w_index, ff_w_index).evaluate(q);
    }
    return ScatteringProfile(Iq, debye_axis);
}

SimpleDataset CompositeDistanceHistogramFFExplicit::debye_transform(const std::vector<double>&) const {
    throw except::not_implemented("CompositeDistanceHistogramFFGrid::debye_transform(const std::vector<double>& q) const");
}

ScatteringProfile CompositeDistanceHistogramFFExplicit::get_profile_ax() const {
    const auto& ff_ax_table = form_factor::storage::cross::get_precalculated_form_factor_table();
    auto sinqd_table = get_sinc_table();
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        double cx = exv_factor(constants::axes::q_vals[q]);
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                double ax_sum = std::inner_product(cp_ax.begin(ff1, ff2), cp_ax.end(ff1, ff2), sinqd_table->begin(q), 0.0);
                Iq[q-q0] += 2*cx*ax_sum*ff_ax_table.index(ff1, ff2).evaluate(q);
            }
        }
    }
    return ScatteringProfile(Iq, debye_axis);
}

ScatteringProfile CompositeDistanceHistogramFFExplicit::get_profile_xx() const {
    const auto& ff_xx_table = form_factor::storage::exv::get_precalculated_form_factor_table();
    auto sinqd_table = get_sinc_table();
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        double cx = exv_factor(constants::axes::q_vals[q]);
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                double xx_sum = std::inner_product(cp_xx.begin(ff1, ff2), cp_xx.end(ff1, ff2), sinqd_table->begin(q), 0.0);
                Iq[q-q0] += cx*cx*xx_sum*ff_xx_table.index(ff1, ff2).evaluate(q);
            }
        }
    }
    return ScatteringProfile(Iq, debye_axis);
}

ScatteringProfile CompositeDistanceHistogramFFExplicit::get_profile_wx() const {
    const auto& ff_ax_table = form_factor::storage::cross::get_precalculated_form_factor_table();
    auto sinqd_table = get_sinc_table();
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    unsigned int ff_w_index = static_cast<int>(form_factor::form_factor_t::OH);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        double cx = exv_factor(constants::axes::q_vals[q]);
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            double wx_sum = std::inner_product(cp_wx.begin(ff1), cp_wx.end(ff1), sinqd_table->begin(q), 0.0);
            Iq[q-q0] += 2*cx*cw*wx_sum*ff_ax_table.index(ff_w_index, ff1).evaluate(q);
        }
    }
    return ScatteringProfile(Iq, debye_axis);
}