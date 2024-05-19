/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hist/intensity_calculator/crysol/CompositeDistanceHistogramCrysol.h>
#include <hist/intensity_calculator/crysol/ExvFormFactorCrysol.h>
#include <settings/HistogramSettings.h>

using namespace hist;

double CompositeDistanceHistogramCrysol::exv_factor(double q) const {
    // G(q) factor from CRYSOL: https://doi.org/10.1107/S0021889895007047
    // double magic_constant = 1/(4*constants::pi*constants::pi);
    double rm = 1.62;
    double c = constexpr_math::pow(4*constants::pi/3, 3./2)*constants::pi*rm*rm;
    return std::pow(cx, 3)*std::exp(-c*(cx*cx - 1)*q*q);
}

Limit CompositeDistanceHistogramCrysol::get_excluded_volume_scaling_factor_limits() const {
    return {0.8, 1.265};
}

SimpleDataset CompositeDistanceHistogramCrysol::debye_transform(const std::vector<double>&) const {
    throw std::runtime_error("CompositeDistanceHistogramCrysol::debye_transform: Custom q values are not supported.");
}

const auto& ff_aa_table = form_factor::storage::atomic::get_precalculated_form_factor_table();
static auto ff_ax_table = form_factor::crysol::storage::cross::generate_table();
static auto ff_xx_table = form_factor::crysol::storage::exv::generate_table();
ScatteringProfile CompositeDistanceHistogramCrysol::debye_transform() const {
    auto sinqd_table = get_sinc_table();
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    unsigned int ff_w_index = static_cast<int>(form_factor::form_factor_t::OH);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        double Gq = exv_factor(constants::axes::q_vals[q]);
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                double count_sum = std::inner_product(cp_xx.begin(ff1, ff2), cp_xx.end(ff1, ff2), sinqd_table->begin(q), 0.0);

                // atom-atom
                Iq[q-q0] += count_sum*ff_aa_table.index(ff1, ff2).evaluate(q);

                // atom-exv
                Iq[q-q0] -= 2*Gq*count_sum*ff_ax_table.index(ff1, ff2).evaluate(q);

                // exv-exv
                Iq[q-q0] += Gq*Gq*count_sum*ff_xx_table.index(ff1, ff2).evaluate(q);
            }

            // the sum is multiplied by the water charge, but this can be absorbed into the cw scaling factor
            double count_sum = std::inner_product(cp_wx.begin(ff1), cp_wx.end(ff1), sinqd_table->begin(q), 0.0);

            // atom-water
            Iq[q-q0] += 2*cw*count_sum*ff_aa_table.index(ff1, ff_w_index).evaluate(q);

            // exv-water
            Iq[q-q0] -= 2*Gq*cw*count_sum*ff_ax_table.index(ff_w_index, ff1).evaluate(q);
        }

        // water-water
        double ww_sum = std::inner_product(cp_ww.begin(), cp_ww.end(), sinqd_table->begin(q), 0.0);
        Iq[q-q0] += cw*cw*ww_sum*ff_aa_table.index(ff_w_index, ff_w_index).evaluate(q);
    }
    return ScatteringProfile(Iq, debye_axis);
}

ScatteringProfile CompositeDistanceHistogramCrysol::get_profile_ax() const {
    auto sinqd_table = get_sinc_table();
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        double Gq = exv_factor(constants::axes::q_vals[q]);
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                double count_sum = std::inner_product(cp_xx.begin(ff1, ff2), cp_xx.end(ff1, ff2), sinqd_table->begin(q), 0.0);
                Iq[q-q0] += 2*Gq*count_sum*ff_ax_table.index(ff1, ff2).evaluate(q);
            }
        }
    }
    return ScatteringProfile(Iq, debye_axis);
}

ScatteringProfile CompositeDistanceHistogramCrysol::get_profile_xx() const {
    auto sinqd_table = get_sinc_table();
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        double Gq = exv_factor(constants::axes::q_vals[q]);
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                double count_sum = std::inner_product(cp_xx.begin(ff1, ff2), cp_xx.end(ff1, ff2), sinqd_table->begin(q), 0.0);
                Iq[q-q0] += Gq*Gq*count_sum*ff_xx_table.index(ff1, ff2).evaluate(q);
            }
        }
    }
    return ScatteringProfile(Iq, debye_axis);}

ScatteringProfile CompositeDistanceHistogramCrysol::get_profile_wx() const {
    auto sinqd_table = get_sinc_table();
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    unsigned int ff_w_index = static_cast<int>(form_factor::form_factor_t::OH);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        double Gq = exv_factor(constants::axes::q_vals[q]);
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            double count_sum = std::inner_product(cp_wx.begin(ff1), cp_wx.end(ff1), sinqd_table->begin(q), 0.0);
            Iq[q-q0] += 2*Gq*cw*count_sum*ff_ax_table.index(ff_w_index, ff1).evaluate(q);
        }
    }
    return ScatteringProfile(Iq, debye_axis);
}

ScatteringProfile CompositeDistanceHistogramCrysol::get_profile_aa() const {
    auto sinqd_table = get_sinc_table();
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                double count_sum = std::inner_product(cp_xx.begin(ff1, ff2), cp_xx.end(ff1, ff2), sinqd_table->begin(q), 0.0);
                Iq[q-q0] += count_sum*ff_aa_table.index(ff1, ff2).evaluate(q);
            }
        }
    }
    return ScatteringProfile(Iq, debye_axis);
}

ScatteringProfile CompositeDistanceHistogramCrysol::get_profile_aw() const {
    auto sinqd_table = get_sinc_table();
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    unsigned int ff_water_index = static_cast<int>(form_factor::form_factor_t::O);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            double count_sum = std::inner_product(cp_wx.begin(ff1), cp_wx.end(ff1), sinqd_table->begin(q), 0.0);
            Iq[q-q0] += 2*cw*count_sum*ff_aa_table.index(ff1, ff_water_index).evaluate(q);
        }
    }
    return ScatteringProfile(Iq, debye_axis);
}

ScatteringProfile CompositeDistanceHistogramCrysol::get_profile_ww() const {
    auto sinqd_table = get_sinc_table();
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    unsigned int ff_water_index = static_cast<int>(form_factor::form_factor_t::O);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        double ww_sum = std::inner_product(cp_ww.begin(), cp_ww.end(), sinqd_table->begin(q), 0.0);
        Iq[q-q0] += cw*cw*ww_sum*ff_aa_table.index(ff_water_index, ff_water_index).evaluate(q);
    }
    return ScatteringProfile(Iq, debye_axis);
}