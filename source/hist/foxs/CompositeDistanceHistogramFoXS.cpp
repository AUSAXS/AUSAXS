#include <hist/foxs/CompositeDistanceHistogramFoXS.h>
#include <hist/foxs/FormFactorFoXS.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFExplicit.h>
#include <hist/Histogram.h>
#include <table/ArrayDebyeTable.h>
#include <settings/HistogramSettings.h>

using namespace hist;

CompositeDistanceHistogramFoXS::CompositeDistanceHistogramFoXS() = default;

CompositeDistanceHistogramFoXS::CompositeDistanceHistogramFoXS(
    container::Container3D<double>&& p_aa, container::Container3D<double>&& p_ax, container::Container3D<double>&& p_xx,
    container::Container2D<double>&& p_wa, container::Container2D<double>&& p_wx, container::Container1D<double>&& p_ww,
    std::vector<double>&& p_tot, const Axis& axis)
: CompositeDistanceHistogramFFAvg(std::move(p_aa), std::move(p_wa), std::move(p_ww), std::move(p_tot), axis), cp_ax(std::move(p_ax)), cp_xx(std::move(p_xx)), cp_wx(std::move(p_wx)) {}

CompositeDistanceHistogramFoXS::~CompositeDistanceHistogramFoXS() = default;

double CompositeDistanceHistogramFoXS::G_factor(double q) const {
    constexpr double rm = 1.58;
    constexpr double c = rm*rm/(4*M_PI);
    // constexpr double c = std::pow(4*M_PI/3, 3./2)*M_PI*1.62*1.62*constants::form_factor::s_to_q_factor;
    return std::pow(cx, 3)*std::exp(-c*(cx*cx - 1)*q*q);
}

static auto ff_aa_table = form_factor::foxs::storage::atomic::generate_table();
static auto ff_ax_table = form_factor::foxs::storage::cross::generate_table();
static auto ff_xx_table = form_factor::foxs::storage::exv::generate_table();
ScatteringProfile CompositeDistanceHistogramFoXS::debye_transform() const {
    const auto& sinqd_table = table::ArrayDebyeTable::get_default_table();

    // calculate the Debye scattering intensity
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    // {
    //     for (unsigned int i = 0; i < 10; ++i) {
    //         double aa_sum = 0;
    //         double xx_sum = 0;
    //         double count = 0;
    //         for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
    //             for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
    //                 count += cp_xx.index(ff1, ff2, i);
    //                 xx_sum += cp_xx.index(ff1, ff2, i)*ff_xx_table.index(ff1, ff2).evaluate(0);
    //                 aa_sum += cp_xx.index(ff1, ff2, i)*ff_aa_table.index(ff1, ff2).evaluate(0);
    //             }
    //         }
    //         std::cout << "d = " << constants::axes::d_axis.get_bin_value(i) << std::endl;
    //         std::cout << "\taa_sum[" << i << "] = " << aa_sum << std::endl;
    //         std::cout << "\txx_sum[" << i << "] = " << xx_sum << std::endl;
    //         std::cout << "\t count[" << i << "] = " << count << std::endl;
    //     }
    // }

    std::vector<double> Iq(debye_axis.bins, 0);
    unsigned int ff_w_index = static_cast<int>(form_factor::form_factor_t::OH);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        double Gq = G_factor(constants::axes::q_vals[q]);
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                double count_sum = std::inner_product(cp_xx.begin(ff1, ff2), cp_xx.end(ff1, ff2), sinqd_table.begin(q), 0.0);

                // atom-atom
                Iq[q] += count_sum*ff_aa_table.index(ff1, ff2).evaluate(q);

                // atom-exv
                Iq[q] -= 2*Gq*count_sum*ff_ax_table.index(ff1, ff2).evaluate(q);

                // exv-exv
                Iq[q] += Gq*Gq*count_sum*ff_xx_table.index(ff1, ff2).evaluate(q);
            }

            // the sum is multiplied by the water charge, but this can be absorbed into the cw scaling factor
            double count_sum = std::inner_product(cp_wx.begin(ff1), cp_wx.end(ff1), sinqd_table.begin(q), 0.0);

            // atom-water
            Iq[q] += 2*cw*count_sum*ff_aa_table.index(ff1, ff_w_index).evaluate(q);

            // exv-water
            Iq[q] -= 2*Gq*cw*count_sum*ff_ax_table.index(ff_w_index, ff1).evaluate(q);
        }

        // water-water
        double ww_sum = std::inner_product(cp_ww.begin(), cp_ww.end(), sinqd_table.begin(q), 0.0);
        Iq[q] += cw*cw*ww_sum*ff_aa_table.index(ff_w_index, ff_w_index).evaluate(q);
    }
    return ScatteringProfile(Iq, debye_axis);
}

const ScatteringProfile CompositeDistanceHistogramFoXS::get_profile_ax() const {
    const auto& sinqd_table = table::ArrayDebyeTable::get_default_table();
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        double Gq = G_factor(constants::axes::q_vals[q]);
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                double count_sum = std::inner_product(cp_xx.begin(ff1, ff2), cp_xx.end(ff1, ff2), sinqd_table.begin(q), 0.0);
                Iq[q] += 2*Gq*count_sum*ff_ax_table.index(ff1, ff2).evaluate(q);
            }
        }
    }
    return ScatteringProfile(Iq, debye_axis);
}

const ScatteringProfile CompositeDistanceHistogramFoXS::get_profile_xx() const {
    const auto& sinqd_table = table::ArrayDebyeTable::get_default_table();
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        double Gq = G_factor(constants::axes::q_vals[q]);
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                double count_sum = std::inner_product(cp_xx.begin(ff1, ff2), cp_xx.end(ff1, ff2), sinqd_table.begin(q), 0.0);
                Iq[q] += Gq*Gq*count_sum*ff_xx_table.index(ff1, ff2).evaluate(q);
            }
        }
    }
    return ScatteringProfile(Iq, debye_axis);}

const ScatteringProfile CompositeDistanceHistogramFoXS::get_profile_wx() const {
    const auto& sinqd_table = table::ArrayDebyeTable::get_default_table();
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    unsigned int ff_w_index = static_cast<int>(form_factor::form_factor_t::OH);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        double Gq = G_factor(constants::axes::q_vals[q]);
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            double count_sum = std::inner_product(cp_wx.begin(ff1), cp_wx.end(ff1), sinqd_table.begin(q), 0.0);
            Iq[q] += 2*Gq*cw*count_sum*ff_ax_table.index(ff_w_index, ff1).evaluate(q);
        }
    }
    return ScatteringProfile(Iq, debye_axis);
}

const ScatteringProfile CompositeDistanceHistogramFoXS::get_profile_aa() const {
    const auto& sinqd_table = table::ArrayDebyeTable::get_default_table();
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                double count_sum = std::inner_product(cp_xx.begin(ff1, ff2), cp_xx.end(ff1, ff2), sinqd_table.begin(q), 0.0);
                Iq[q] += count_sum*ff_aa_table.index(ff1, ff2).evaluate(q);
            }
        }
    }
    return ScatteringProfile(Iq, debye_axis);
}

const ScatteringProfile CompositeDistanceHistogramFoXS::get_profile_aw() const {
    const auto& sinqd_table = table::ArrayDebyeTable::get_default_table();
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    unsigned int ff_water_index = static_cast<int>(form_factor::form_factor_t::O);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            double count_sum = std::inner_product(cp_wx.begin(ff1), cp_wx.end(ff1), sinqd_table.begin(q), 0.0);
            Iq[q] += 2*cw*count_sum*ff_aa_table.index(ff1, ff_water_index).evaluate(q);
        }
    }
    return ScatteringProfile(Iq, debye_axis);
}

const ScatteringProfile CompositeDistanceHistogramFoXS::get_profile_ww() const {
    const auto& sinqd_table = table::ArrayDebyeTable::get_default_table();
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    unsigned int ff_water_index = static_cast<int>(form_factor::form_factor_t::O);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        double ww_sum = std::inner_product(cp_ww.begin(), cp_ww.end(), sinqd_table.begin(q), 0.0);
        Iq[q] += cw*cw*ww_sum*ff_aa_table.index(ff_water_index, ff_water_index).evaluate(q);
    }
    return ScatteringProfile(Iq, debye_axis);
}