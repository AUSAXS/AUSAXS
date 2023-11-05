#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvgBase.h>
#include <hist/Histogram.h>
#include <table/ArrayDebyeTable.h>
#include <form_factor/FormFactor.h>
#include <form_factor/PrecalculatedFormFactorProduct.h>
#include <settings/HistogramSettings.h>

using namespace hist;

template<typename FormFactorTableType>
CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::CompositeDistanceHistogramFFAvgBase() = default;

template<typename FormFactorTableType>
CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::CompositeDistanceHistogramFFAvgBase(
    hist::WeightedDistribution3D&& p_aa, 
    hist::WeightedDistribution2D&& p_aw, 
    hist::WeightedDistribution1D&& p_ww, 
    const Axis& axis
) : ICompositeDistanceHistogramExv(hist::Distribution1D(axis.bins, 0), axis), cp_aa(std::move(p_aa)), cp_aw(std::move(p_aw)), cp_ww(std::move(p_ww)) {
    use_weighted_sinc_table();
}

template<typename FormFactorTableType>
CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::CompositeDistanceHistogramFFAvgBase(
    hist::Distribution3D&& p_aa, 
    hist::Distribution2D&& p_aw, 
    hist::Distribution1D&& p_ww, 
    const Axis& axis
) : ICompositeDistanceHistogramExv(hist::Distribution1D(axis.bins, 0), axis), cp_aa(std::move(p_aa)), cp_aw(std::move(p_aw)), cp_ww(std::move(p_ww)) {}

template<typename FormFactorTableType>
CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::~CompositeDistanceHistogramFFAvgBase() = default;

// #define DEBUG_DEBYE_TRANSFORM 0
// #define DEBUG_PLOT_FF 0
// #if DEBUG_DEBYE_TRANSFORM
//     #include <plots/PlotDataset.h>
//     #include <settings/GeneralSettings.h>
//     #include <dataset/SimpleDataset.h>
//     static unsigned int qcheck = 0;
// #endif
// template<typename FormFactorTableType>
// ScatteringProfile CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::debye_transform() const {
//     const auto& sinqd_table = get_sinc_table();
//     const auto& ff_table = get_ff_table();

//     // calculate the Debye scattering intensity
//     Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
//     unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

//     std::vector<double> Iq(debye_axis.bins, 0);
//     std::vector<double> q_axis = debye_axis.as_vector();

//     #ifdef DEBUG_PLOT_FF
//         std::vector<double> I_aa(q_axis.size(), 0);
//         std::vector<double> I_aw(q_axis.size(), 0);
//         std::vector<double> I_ax(q_axis.size(), 0);
//         std::vector<double> I_ww(q_axis.size(), 0);
//         std::vector<double> I_wx(q_axis.size(), 0);
//         std::vector<double> I_xx(q_axis.size(), 0);
//     #endif

//     for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
//         for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
//             // atom-atom
//             for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
//                 double atom_atom_sum = std::inner_product(cp_aa.begin(ff1, ff2), cp_aa.end(ff1, ff2), sinqd_table.begin(q), 0.0);
//                 Iq[q] += atom_atom_sum*ff_table.index(ff1, ff2).evaluate(q);
//                 #if DEBUG_DEBYE_TRANSFORM
//                     if (q==qcheck && atom_atom_sum != 0) {
//                         std::cout << "(aa) Iq[" << q << "] += aa_sum*ff_table[" << ff1 << ", " << ff2 << "](" << q << ") = " 
//                             << atom_atom_sum << "*" << ff_table.index(ff1, ff2).evaluate(q) << " = " 
//                             << atom_atom_sum*ff_table.index(ff1, ff2).evaluate(q) << std::endl;
//                         for (unsigned int d = 0; d < axis.bins; ++d) {
//                             if (cp_aa.index(ff1, ff2, d) != 0) {std::cout << "\t\taa_sum += p_pp[" << ff1 << ", " << ff2 << ", " << d << "] = " << cp_aa.index(ff1, ff2, d) << "*" << sinqd_table.lookup(q, d) << std::endl;}
//                         }
//                         std::cout << "\taa_sum = " << atom_atom_sum << std::endl;
//                     }
//                     if (q == qcheck && atom_atom_sum != 0) {std::cout << "\tIq[" << q << "] = " << Iq[q] << std::endl;}
//                 #endif
//                 #if DEBUG_PLOT_FF
//                     I_aa[q] += atom_atom_sum*ff_table.index(ff1, ff2).evaluate(q);
//                 #endif
//             }

//             // atom-exv
//             double atom_exv_sum = std::inner_product(cp_aa.begin(ff1, form_factor::exv_bin), cp_aa.end(ff1, form_factor::exv_bin), sinqd_table.begin(q), 0.0);
//             Iq[q] -= 2*cx*atom_exv_sum*ff_table.index(ff1, form_factor::exv_bin).evaluate(q);
//             #if DEBUG_DEBYE_TRANSFORM
//                 if (q==qcheck && atom_exv_sum != 0) {
//                     std::cout << "(ax) Iq[" << q << "] -= 2*ax_sum*ff_table[" << ff1 << ", " << form_factor::exv_bin << "](" << q << ") = " 
//                         << "2*" << atom_exv_sum << "*" << ff_table.index(ff1, form_factor::exv_bin).evaluate(q) << " = " 
//                         << 2*atom_exv_sum*ff_table.index(ff1, form_factor::exv_bin).evaluate(q) << std::endl;
//                     for (unsigned int d = 0; d < axis.bins; ++d) {
//                         if (cp_aa.index(ff1, form_factor::exv_bin, d) != 0) {std::cout << "\t\tax_sum += 2*p_pp[" << ff1 << ", " << form_factor::exv_bin << ", " << d << "] = " << 2*(cp_aa.index(ff1, form_factor::exv_bin, d) + cp_aa.index(form_factor::exv_bin, ff1, d)) << "*" << sinqd_table.lookup(q, d) << std::endl;}
//                     }
//                     std::cout << "\tax_sum = " << atom_exv_sum << std::endl;
//                 }
//                 if (q == qcheck && atom_exv_sum != 0) {std::cout << "\tIq[" << q << "] = " << Iq[q] << std::endl;}
//             #endif
//             #if DEBUG_PLOT_FF
//                 I_ax[q] -= 2*cx*atom_exv_sum*ff_table.index(ff1, form_factor::exv_bin).evaluate(q);
//             #endif

//             // atom-water
//             double atom_water_sum = std::inner_product(cp_aw.begin(ff1), cp_aw.end(ff1), sinqd_table.begin(q), 0.0);
//             Iq[q] += 2*cw*atom_water_sum*ff_table.index(ff1, form_factor::water_bin).evaluate(q);
//             #if DEBUG_DEBYE_TRANSFORM
//                 if (q==qcheck && ff1 == 1) {
//                     std::cout << "(aw) Iq[" << q << "] += 2*aw_sum*ff_table[" << ff1 << ", " << form_factor::water_bin << "](" << q << ") = " 
//                         << "2*" << atom_water_sum << "*" << ff_table.index(ff1, form_factor::water_bin).evaluate(q) << " = "
//                         << 2*atom_water_sum*ff_table.index(ff1, form_factor::water_bin).evaluate(q) << std::endl;
//                     for (unsigned int d = 0; d < axis.bins; ++d) {
//                         if (cp_aw.index(ff1, d) != 0) {std::cout << "\t\taw_sum += p_hp[" << ff1 << ", " << d << "] = " << cp_aw.index(ff1, d) << "*" << sinqd_table.lookup(q, d) << " = " << cp_aw.index(ff1, d)*sinqd_table.lookup(q, d) << std::endl;}
//                     }
//                     std::cout << "\taw_sum = " << atom_water_sum << std::endl;
//                 }
//                 if (q == qcheck && ff1 == 1) {std::cout << "\tIq[" << q << "] = " << Iq[q] << std::endl;}
//             #endif
//             #if DEBUG_PLOT_FF
//                 I_aw[q] += 2*cw*atom_water_sum*ff_table.index(ff1, form_factor::water_bin).evaluate(q);
//             #endif
//         }

//         // exv-exv
//         double exv_exv_sum = std::inner_product(cp_aa.begin(form_factor::exv_bin, form_factor::exv_bin), cp_aa.end(form_factor::exv_bin, form_factor::exv_bin), sinqd_table.begin(q), 0.0);
//         Iq[q] += cx*cx*exv_exv_sum*ff_table.index(form_factor::exv_bin, form_factor::exv_bin).evaluate(q);
//         #if DEBUG_DEBYE_TRANSFORM
//             if (q==qcheck) {
//                 std::cout << "(xx) Iq[" << q << "] += xx_sum*ff_table[" << form_factor::exv_bin << ", " << form_factor::exv_bin << "](" << q << ") = " 
//                     << exv_exv_sum << "*" << ff_table.index(form_factor::exv_bin, form_factor::exv_bin).evaluate(q) << " = " 
//                     << exv_exv_sum*ff_table.index(form_factor::exv_bin, form_factor::exv_bin).evaluate(q) << std::endl;
//                 for (unsigned int d = 0; d < axis.bins; ++d) {
//                     if (cp_aa.index(form_factor::exv_bin, form_factor::exv_bin, d) != 0) {std::cout << "\t\txx_sum += p_pp[" << form_factor::exv_bin << ", " << form_factor::exv_bin << ", " << d << "] = " << cp_aa.index(form_factor::exv_bin, form_factor::exv_bin, d) << std::endl;}
//                 }
//                 std::cout << "\txx_sum = " << exv_exv_sum << std::endl;
//             }
//             if (q == qcheck) {std::cout << "\tIq[" << q << "] = " << Iq[q] << std::endl;}
//         #endif
//         #if DEBUG_PLOT_FF
//             I_xx[q] += cx*cx*exv_exv_sum*ff_table.index(form_factor::exv_bin, form_factor::exv_bin).evaluate(q);
//         #endif

//         // exv-water
//         double exv_water_sum = std::inner_product(cp_aw.begin(form_factor::exv_bin), cp_aw.end(form_factor::exv_bin), sinqd_table.begin(q), 0.0);
//         Iq[q] -= 2*cx*cw*exv_water_sum*ff_table.index(form_factor::exv_bin, form_factor::water_bin).evaluate(q);
//         #if DEBUG_DEBYE_TRANSFORM
//             if (q==qcheck) {
//                 std::cout << "(wx) Iq[" << q << "] -= 2*wx_sum*ff_table[" << form_factor::exv_bin << ", " << form_factor::water_bin << "](" << q << ") = "
//                      << "2*" << exv_water_sum << "*" << ff_table.index(form_factor::exv_bin, form_factor::water_bin).evaluate(q) << " = " 
//                      << 2*exv_water_sum*ff_table.index(form_factor::exv_bin, form_factor::water_bin).evaluate(q) << std::endl;
//                 for (unsigned int d = 0; d < axis.bins; ++d) {
//                     if (cp_aw.index(form_factor::exv_bin, d) != 0) {std::cout << "\t\twx_sum += p_hp[" << form_factor::exv_bin << ", " << d << "] = " << cp_aw.index(form_factor::exv_bin, d) << "*" << sinqd_table.lookup(q, d) << std::endl;}
//                 }
//                 std::cout << "\twx_sum = " << exv_water_sum << std::endl;
//             }
//             if (q == qcheck) {std::cout << "\tIq[" << q << "] = " << Iq[q] << std::endl;}
//         #endif
//         #if DEBUG_PLOT_FF
//             I_wx[q] -= 2*cx*cw*exv_water_sum*ff_table.index(form_factor::exv_bin, form_factor::water_bin).evaluate(q);
//         #endif

//         // water-water
//         double water_water_sum = std::inner_product(cp_ww.begin(), cp_ww.end(), sinqd_table.begin(q), 0.0);
//         Iq[q] += cw*cw*water_water_sum*ff_table.index(form_factor::water_bin, form_factor::water_bin).evaluate(q);
//         #if DEBUG_DEBYE_TRANSFORM
//             if (q==qcheck) {std::cout << "(ww) Iq[" << q << "] += ww_sum*ff_table[" << form_factor::water_bin << ", " << form_factor::water_bin << "](" << q << ") = " 
//                 << water_water_sum << "*" << ff_table.index(form_factor::water_bin, form_factor::water_bin).evaluate(q) << " = " 
//                 << water_water_sum*ff_table.index(form_factor::water_bin, form_factor::water_bin).evaluate(q) << std::endl;}
//             if (q == qcheck) {std::cout << "\tIq[" << q << "] = " << Iq[q] << std::endl;}
//         #endif
//         #if DEBUG_PLOT_FF
//             I_ww[q] += cw*cw*water_water_sum*ff_table.index(form_factor::water_bin, form_factor::water_bin).evaluate(q);
//         #endif
//     }

//     #if DEBUG_PLOT_FF
//         for (unsigned int i = 0; i < q_axis.size(); ++i) {
//             I_aa[i] = std::abs(I_aa[i]);
//             I_aw[i] = std::abs(I_aw[i]);
//             I_ax[i] = std::abs(I_ax[i]);
//             I_ww[i] = std::abs(I_ww[i]);
//             I_wx[i] = std::abs(I_wx[i]);
//             I_xx[i] = std::abs(I_xx[i]);
//         }
//         SimpleDataset temp_aa(q_axis, I_aa);
//         SimpleDataset temp_aw(q_axis, I_aw);
//         SimpleDataset temp_ax(q_axis, I_ax);
//         SimpleDataset temp_ww(q_axis, I_ww);
//         SimpleDataset temp_wx(q_axis, I_wx);
//         SimpleDataset temp_xx(q_axis, I_xx);
//         temp_aa.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"ylabel", "I"}, {"legend", "I_aa"}, {"color", style::color::red}});
//         temp_aw.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"ylabel", "I"}, {"legend", "I_aw"}, {"color", style::color::green}});
//         temp_ax.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"ylabel", "I"}, {"legend", "I_ax"}, {"color", style::color::blue}});
//         temp_ww.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"ylabel", "I"}, {"legend", "I_ww"}, {"color", style::color::brown}});
//         temp_wx.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"ylabel", "I"}, {"legend", "I_wx"}, {"color", style::color::orange}});
//         temp_xx.add_plot_options({{"xlabel", "q"}, {"linewidth", 2}, {"logx", true}, {"logy", true}, {"ylabel", "I"}, {"legend", "I_xx"}, {"color", style::color::purple}});

//         plots::PlotDataset(temp_aa)
//             .plot(temp_ax)
//             .plot(temp_xx)
//         .save(settings::general::output + "ff.png");
//         temp_aa.save(settings::general::output + "ff_aa.dat");
//         temp_ax.save(settings::general::output + "ff_ax.dat");
//         temp_xx.save(settings::general::output + "ff_xx.dat");
//         temp_aw.save(settings::general::output + "ff_aw.dat");
//         temp_wx.save(settings::general::output + "ff_wx.dat");
//         temp_ww.save(settings::general::output + "ff_ww.dat");
//     #endif

//     return ScatteringProfile(Iq, debye_axis);
// }

template<typename FormFactorTableType>
ScatteringProfile CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::debye_transform() const {
    const auto& ff_table = get_ff_table();
    auto sinqd_table = get_sinc_table();

    // calculate the Debye scattering intensity
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            // atom-atom
            for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                double aa_sum = std::inner_product(cp_aa.begin(ff1, ff2), cp_aa.end(ff1, ff2), sinqd_table->begin(q), 0.0);
                Iq[q] += aa_sum*ff_table.index(ff1, ff2).evaluate(q);
            }

            // atom-exv
            double ax_sum = std::inner_product(cp_aa.begin(ff1, form_factor::exv_bin), cp_aa.end(ff1, form_factor::exv_bin), sinqd_table->begin(q), 0.0);
            Iq[q] -= 2*cx*ax_sum*ff_table.index(ff1, form_factor::exv_bin).evaluate(q);

            // atom-water
            double aw_sum = std::inner_product(cp_aw.begin(ff1), cp_aw.end(ff1), sinqd_table->begin(q), 0.0);
            Iq[q] += 2*cw*aw_sum*ff_table.index(ff1, form_factor::water_bin).evaluate(q);
        }

        // exv-exv
        double xx_sum = std::inner_product(cp_aa.begin(form_factor::exv_bin, form_factor::exv_bin), cp_aa.end(form_factor::exv_bin, form_factor::exv_bin), sinqd_table->begin(q), 0.0);
        Iq[q] += cx*cx*xx_sum*ff_table.index(form_factor::exv_bin, form_factor::exv_bin).evaluate(q);

        // exv-water
        double ew_sum = std::inner_product(cp_aw.begin(form_factor::exv_bin), cp_aw.end(form_factor::exv_bin), sinqd_table->begin(q), 0.0);
        Iq[q] -= 2*cx*cw*ew_sum*ff_table.index(form_factor::exv_bin, form_factor::water_bin).evaluate(q);

        // water-water
        double ww_sum = std::inner_product(cp_ww.begin(), cp_ww.end(), sinqd_table->begin(q), 0.0);
        Iq[q] += cw*cw*ww_sum*ff_table.index(form_factor::water_bin, form_factor::water_bin).evaluate(q);
    }
    return ScatteringProfile(Iq, debye_axis);
}

template<typename FormFactorTableType>
const std::vector<double>& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_counts() const {
    p = std::vector<double>(axis.bins, 0);
    auto& p_pp = get_aa_counts();
    auto& p_hp = get_aw_counts();
    auto& p_hh = get_ww_counts();
    for (unsigned int i = 0; i < axis.bins; ++i) {
        p[i] = p_pp.index(i) + 2*cw*p_hp.index(i) + cw*cw*p_hh.index(i);
    }
    return p.data;
}

template<typename FormFactorTableType>
Distribution1D& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_aa_counts() {
    return const_cast<Distribution1D&>(const_cast<const CompositeDistanceHistogramFFAvgBase*>(this)->get_aa_counts());
}

template<typename FormFactorTableType>
const Distribution1D& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_aa_counts() const {
    if (!p_aa.empty()) {return p_aa;}
    p_aa = Distribution1D(axis.bins, 0);
    for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
        for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
            std::transform(p_aa.begin(), p_aa.end(), cp_aa.begin(ff1, ff2), p_aa.begin(), std::plus<>());
        }
    }
    return p_aa;
}

template<typename FormFactorTableType>
Distribution1D& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_aw_counts() {return const_cast<Distribution1D&>(const_cast<const CompositeDistanceHistogramFFAvgBase*>(this)->get_aw_counts());}

template<typename FormFactorTableType>
const Distribution1D& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_aw_counts() const {
    if (!p_aw.empty()) {return p_aw;}
    p_aw = Distribution1D(axis.bins, 0);
    for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
        std::transform(p_aw.begin(), p_aw.end(), cp_aw.begin(ff1), p_aw.begin(), std::plus<>());
    }
    return p_aw;
}

template<typename FormFactorTableType>
Distribution1D& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_ww_counts() {return const_cast<Distribution1D&>(const_cast<const CompositeDistanceHistogramFFAvgBase*>(this)->get_ww_counts());}

template<typename FormFactorTableType>
const Distribution1D& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_ww_counts() const {
    if (!p_ww.empty()) {return p_ww;}
    p_ww = Distribution1D(axis.bins, 0);
    std::transform(p_ww.begin(), p_ww.end(), cp_ww.begin(), p_ww.begin(), std::plus<>());
    return p_ww;
}

template<typename FormFactorTableType>
const Distribution3D& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_aa_counts_ff() const {
    return cp_aa;
}

template<typename FormFactorTableType>
Distribution3D& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_aa_counts_ff() {
    return cp_aa;
}

template<typename FormFactorTableType>
const Distribution2D& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_aw_counts_ff() const {
    return cp_aw;
}

template<typename FormFactorTableType>
Distribution2D& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_aw_counts_ff() {
    return cp_aw;
}

template<typename FormFactorTableType>
const Distribution1D& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_ww_counts_ff() const {
    return cp_ww;
}

template<typename FormFactorTableType>
Distribution1D& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_ww_counts_ff() {
    return cp_ww;
}

template<typename FormFactorTableType>
void CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::apply_water_scaling_factor(double k) {
    cw = k;
}

template<typename FormFactorTableType>
void CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::apply_excluded_volume_scaling_factor(double k) {
    cx = k;
}

template<typename FormFactorTableType>
const ScatteringProfile CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_profile_aa() const {
    const auto& ff_table = get_ff_table();
    auto sinqd_table = get_sinc_table();
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                double aa_sum = std::inner_product(cp_aa.begin(ff1, ff2), cp_aa.end(ff1, ff2), sinqd_table->begin(q), 0.0);
                Iq[q] += aa_sum*ff_table.index(ff1, ff2).evaluate(q);
            }
        }
    }
    return ScatteringProfile(Iq, debye_axis);
}

template<typename FormFactorTableType>
const ScatteringProfile CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_profile_ax() const {
    const auto& ff_table = get_ff_table();
    auto sinqd_table = get_sinc_table();
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            double ax_sum = std::inner_product(cp_aa.begin(ff1, form_factor::exv_bin), cp_aa.end(ff1, form_factor::exv_bin), sinqd_table->begin(q), 0.0);
            Iq[q] += 2*cx*ax_sum*ff_table.index(ff1, form_factor::exv_bin).evaluate(q);
        }
    }
    return ScatteringProfile(Iq, debye_axis);
}

template<typename FormFactorTableType>
const ScatteringProfile CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_profile_xx() const {
    const auto& ff_table = get_ff_table();
    auto sinqd_table = get_sinc_table();
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        double xx_sum = std::inner_product(cp_aa.begin(form_factor::exv_bin, form_factor::exv_bin), cp_aa.end(form_factor::exv_bin, form_factor::exv_bin), sinqd_table->begin(q), 0.0);
        Iq[q] += cx*cx*xx_sum*ff_table.index(form_factor::exv_bin, form_factor::exv_bin).evaluate(q);
    }
    return ScatteringProfile(Iq, debye_axis);
}

template<typename FormFactorTableType>
const ScatteringProfile CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_profile_wx() const {
    const auto& ff_table = get_ff_table();
    auto sinqd_table = get_sinc_table();
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        double ew_sum = std::inner_product(cp_aw.begin(form_factor::exv_bin), cp_aw.end(form_factor::exv_bin), sinqd_table->begin(q), 0.0);
        Iq[q] += 2*cx*cw*ew_sum*ff_table.index(form_factor::exv_bin, form_factor::water_bin).evaluate(q);
    }
    return ScatteringProfile(Iq, debye_axis);
}

template<typename FormFactorTableType>
const ScatteringProfile CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_profile_aw() const {
    const auto& ff_table = get_ff_table();
    auto sinqd_table = get_sinc_table();
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            double aw_sum = std::inner_product(cp_aw.begin(ff1), cp_aw.end(ff1), sinqd_table->begin(q), 0.0);
            Iq[q] += 2*cw*aw_sum*ff_table.index(ff1, form_factor::water_bin).evaluate(q);
        }
    }
    return ScatteringProfile(Iq, debye_axis);
}

template<typename FormFactorTableType>
const ScatteringProfile CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_profile_ww() const {
    const auto& ff_table = get_ff_table();
    auto sinqd_table = get_sinc_table();
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        double ww_sum = std::inner_product(cp_ww.begin(), cp_ww.end(), sinqd_table->begin(q), 0.0);
        Iq[q] += cw*cw*ww_sum*ff_table.index(form_factor::water_bin, form_factor::water_bin).evaluate(q);
    }
    return ScatteringProfile(Iq, debye_axis);
}

template class hist::CompositeDistanceHistogramFFAvgBase<form_factor::storage::atomic::table_t>;