#include <hist/CompositeDistanceHistogramFF.h>
#include <hist/CompositeDistanceHistogram.h>
#include <hist/Histogram.h>
#include <table/ArrayDebyeTable.h>
#include <form_factor/FormFactor.h>
#include <form_factor/PrecalculatedFormFactorProduct.h>
#include <dataset/SimpleDataset.h>
#include <data/Molecule.h>
#include <data/record/Atom.h>
#include <settings/HistogramSettings.h>

using namespace hist;

CompositeDistanceHistogramFF::CompositeDistanceHistogramFF(container::Container3D<double>&& p_pp, container::Container2D<double>&& p_hp, container::Container1D<double>&& p_hh, std::vector<double>&& p_tot, const Axis& axis) 
    : CompositeDistanceHistogram(std::move(p_tot), axis), cp_pp(std::move(p_pp)), cp_hp(std::move(p_hp)), cp_hh(std::move(p_hh)) {}

CompositeDistanceHistogramFF::CompositeDistanceHistogramFF(CompositeDistanceHistogramFF&& other) noexcept 
    : CompositeDistanceHistogram(std::move(other.p_pp), std::move(other.p_hp), std::move(other.p_hh), std::move(other.p), other.axis), w_scaling(other.w_scaling), exv_scaling(other.exv_scaling), cp_pp(std::move(other.cp_pp)), cp_hp(std::move(other.cp_hp)), cp_hh(std::move(other.cp_hh)) {}

CompositeDistanceHistogramFF::~CompositeDistanceHistogramFF() = default;

#include <plots/PlotDataset.h>
#include <settings/GeneralSettings.h>

#define DEBUG_DEBYE_TRANSFORM 0
#define DEBUG_PLOT_FF 0
#if DEBUG_DEBYE_TRANSFORM
    static unsigned int qcheck = 1;
#endif

// ScatteringProfile CompositeDistanceHistogramFF::debye_transform() const {
//     static Container2D<form_factor::PrecalculatedFormFactorProduct> ff_table = form_factor::PrecalculatedFormFactorProduct::generate_table();

//     // calculate the Debye scattering intensity
//     Axis debye_axis(settings::axes::qmin, settings::axes::qmax, settings::axes::bins);

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

//     unsigned int excluded_volume_index = static_cast<int>(form_factor::form_factor_t::EXCLUDED_VOLUME);
//     unsigned int neutral_oxygen_index = static_cast<int>(form_factor::form_factor_t::O);

//     for (unsigned int q = 0; q < debye_axis.bins; ++q) {
//         for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
//             // atom-atom
//             for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
//                 double atom_atom_sum = 0;
//                 for (unsigned int d = 0; d < axis.bins; ++d) {
//                     atom_atom_sum += cp_pp.index(ff1, ff2, d)*sinqd_table->lookup(q, d);
//                 }
//                 Iq[q] += atom_atom_sum*ff_table.index(ff1, ff2).evaluate(q);
//                 #if DEBUG_DEBYE_TRANSFORM
//                     if (q==qcheck && atom_atom_sum != 0) {
//                         std::cout << "(aa) Iq[" << q << "] += aa_sum*ff_table[" << ff1 << ", " << ff2 << "](" << q << ") = " 
//                             << atom_atom_sum << "*" << ff_table.index(ff1, ff2).evaluate(q) << " = " 
//                             << atom_atom_sum*ff_table.index(ff1, ff2).evaluate(q) << std::endl;
//                         for (unsigned int d = 0; d < axis.bins; ++d) {
//                             if (cp_pp.index(ff1, ff2, d) != 0) {std::cout << "\t\taa_sum += p_pp[" << ff1 << ", " << ff2 << ", " << d << "] = " << cp_pp.index(ff1, ff2, d) << "*" << sinqd_table->lookup(q, d) << std::endl;}
//                         }
//                         std::cout << "\taa_sum = " << atom_atom_sum << std::endl;
//                     }
//                     if (q == qcheck && atom_atom_sum != 0) {std::cout << "\tIq[" << q << "] = " << Iq[q] << std::endl;}
//                 #endif
//                 #if DEBUG_PLOT_FF
//                     I_aa[q] += atom_atom_sum*ff_table.index(ff1, ff2).evaluate(q);
//                 #endif
//             }

//             double atom_exv_sum = 0; // atom-exv
//             for (unsigned int d = 0; d < axis.bins; ++d) {
//                 atom_exv_sum += cp_pp.index(ff1, excluded_volume_index, d)*sinqd_table->lookup(q, d);
//             }
//             Iq[q] -= 2*exv_scaling*atom_exv_sum*ff_table.index(ff1, excluded_volume_index).evaluate(q);
//             #if DEBUG_DEBYE_TRANSFORM
//                 if (q==qcheck && atom_exv_sum != 0) {
//                     std::cout << "(ae) Iq[" << q << "] -= 2*ax_sum*ff_table[" << ff1 << ", " << excluded_volume_index << "](" << q << ") = " 
//                         << "2*" << atom_exv_sum << "*" << ff_table.index(ff1, excluded_volume_index).evaluate(q) << " = " 
//                         << 2*atom_exv_sum*ff_table.index(ff1, excluded_volume_index).evaluate(q) << std::endl;
//                     for (unsigned int d = 0; d < axis.bins; ++d) {
//                         if (cp_pp.index(ff1, excluded_volume_index, d) != 0) {std::cout << "\t\tax_sum += 2*p_pp[" << ff1 << ", " << excluded_volume_index << ", " << d << "] = " << 2*(cp_pp.index(ff1, excluded_volume_index, d) + cp_pp.index(excluded_volume_index, ff1, d)) << "*" << sinqd_table->lookup(q, d) << std::endl;}
//                     }
//                     std::cout << "\tax_sum = " << atom_exv_sum << std::endl;
//                 }
//                 if (q == qcheck && atom_exv_sum != 0) {std::cout << "\tIq[" << q << "] = " << Iq[q] << std::endl;}
//             #endif
//             #if DEBUG_PLOT_FF
//                 I_ax[q] -= 2*exv_scaling*atom_exv_sum*ff_table.index(ff1, excluded_volume_index).evaluate(q);
//             #endif

//             // atom-water
//             double atom_water_sum = 0;
//             for (unsigned int d = 0; d < axis.bins; ++d) {
//                 atom_water_sum += cp_hp.index(ff1, d)*sinqd_table->lookup(q, d);
//             }
//             Iq[q] += 2*w_scaling*atom_water_sum*ff_table.index(ff1, neutral_oxygen_index).evaluate(q);
//             #if DEBUG_DEBYE_TRANSFORM
//                 if (q==qcheck && ff1 == 1) {
//                     std::cout << "(aw) Iq[" << q << "] += 2*aw_sum*ff_table[" << ff1 << ", " << neutral_oxygen_index << "](" << q << ") = " 
//                         << "2*" << atom_water_sum << "*" << ff_table.index(ff1, neutral_oxygen_index).evaluate(q) 
//                         << 2*atom_water_sum*ff_table.index(ff1, neutral_oxygen_index).evaluate(q) << std::endl;
//                     for (unsigned int d = 0; d < axis.bins; ++d) {
//                         if (cp_hp.index(ff1, d) != 0) {std::cout << "\t\taw_sum += p_hp[" << ff1 << ", " << d << "] = " << cp_hp.index(ff1, d) << "*" << sinqd_table->lookup(q, d) << " = " << cp_hp.index(ff1, d)*sinqd_table->lookup(q, d) << std::endl;}
//                     }
//                     std::cout << "\taw_sum = " << atom_water_sum << std::endl;
//                 }
//                 if (q == qcheck && ff1 == 1) {std::cout << "\tIq[" << q << "] = " << Iq[q] << std::endl;}
//             #endif
//             #if DEBUG_PLOT_FF
//                 I_aw[q] += 2*w_scaling*atom_water_sum*ff_table.index(ff1, neutral_oxygen_index).evaluate(q);
//             #endif
//         }

//         // exv-exv
//         double exv_exv_sum = 0;
//         for (unsigned int d = 0; d < axis.bins; ++d) {
//             exv_exv_sum += cp_pp.index(excluded_volume_index, excluded_volume_index, d)*sinqd_table->lookup(q, d);
//         }
//         Iq[q] += exv_scaling*exv_scaling*exv_exv_sum*ff_table.index(excluded_volume_index, excluded_volume_index).evaluate(q);
//         #if DEBUG_DEBYE_TRANSFORM
//             if (q==qcheck) {
//                 std::cout << "(ee) Iq[" << q << "] += xx_sum*ff_table[" << excluded_volume_index << ", " << excluded_volume_index << "](" << q << ") = " 
//                     << exv_exv_sum << "*" << ff_table.index(excluded_volume_index, excluded_volume_index).evaluate(q) << " = " 
//                     << exv_exv_sum*ff_table.index(excluded_volume_index, excluded_volume_index).evaluate(q) << std::endl;
//                 for (unsigned int d = 0; d < axis.bins; ++d) {
//                     if (cp_pp.index(excluded_volume_index, excluded_volume_index, d) != 0) {std::cout << "\t\txx_sum += p_pp[" << excluded_volume_index << ", " << excluded_volume_index << ", " << d << "] = " << cp_pp.index(excluded_volume_index, excluded_volume_index, d) << std::endl;}
//                 }
//                 std::cout << "\txx_sum = " << exv_exv_sum << std::endl;
//             }
//             if (q == qcheck) {std::cout << "\tIq[" << q << "] = " << Iq[q] << std::endl;}
//         #endif
//         #if DEBUG_PLOT_FF
//             I_xx[q] += exv_scaling*exv_scaling*exv_exv_sum*ff_table.index(excluded_volume_index, excluded_volume_index).evaluate(q);
//         #endif

//         // exv-water
//         double exv_water_sum = 0;
//         for (unsigned int d = 0; d < axis.bins; ++d) {
//             exv_water_sum += cp_hp.index(excluded_volume_index, d)*sinqd_table->lookup(q, d);
//         }
//         Iq[q] -= 2*exv_scaling*w_scaling*exv_water_sum*ff_table.index(excluded_volume_index, neutral_oxygen_index).evaluate(q);
//         #if DEBUG_DEBYE_TRANSFORM
//             if (q==qcheck) {
//                 std::cout << "(ew) Iq[" << q << "] -= 2*wx_sum*ff_table[" << excluded_volume_index << ", " << neutral_oxygen_index << "](" << q << ") = "
//                      << "2*" << exv_water_sum << "*" << ff_table.index(excluded_volume_index, neutral_oxygen_index).evaluate(q) << " = " 
//                      << 2*exv_water_sum*ff_table.index(excluded_volume_index, neutral_oxygen_index).evaluate(q) << std::endl;
//                 for (unsigned int d = 0; d < axis.bins; ++d) {
//                     if (cp_hp.index(excluded_volume_index, d) != 0) {std::cout << "\t\twx_sum += p_hp[" << excluded_volume_index << ", " << d << "] = " << cp_hp.index(excluded_volume_index, d) << "*" << sinqd_table->lookup(q, d) << std::endl;}
//                 }
//                 std::cout << "\twx_sum = " << exv_water_sum << std::endl;
//             }
//             if (q == qcheck) {std::cout << "\tIq[" << q << "] = " << Iq[q] << std::endl;}
//         #endif
//         #if DEBUG_PLOT_FF
//             I_wx[q] -= 2*exv_scaling*w_scaling*exv_water_sum*ff_table.index(excluded_volume_index, neutral_oxygen_index).evaluate(q);
//         #endif

//         // water-water
//         double water_water_sum = 0;
//         for (unsigned int d = 0; d < axis.bins; ++d) {
//             water_water_sum += cp_hh.index(d)*sinqd_table->lookup(q, d);
//         }
//         Iq[q] += w_scaling*w_scaling*water_water_sum*ff_table.index(neutral_oxygen_index, neutral_oxygen_index).evaluate(q);
//         #if DEBUG_DEBYE_TRANSFORM
//             if (q==qcheck) {std::cout << "(ww) Iq[" << q << "] += ww_sum*ff_table[" << neutral_oxygen_index << ", " << neutral_oxygen_index << "](" << q << ") = " 
//                 << water_water_sum << "*" << ff_table.index(neutral_oxygen_index, neutral_oxygen_index).evaluate(q) << " = " 
//                 << water_water_sum*ff_table.index(neutral_oxygen_index, neutral_oxygen_index).evaluate(q) << std::endl;}
//             if (q == qcheck) {std::cout << "\tIq[" << q << "] = " << Iq[q] << std::endl;}
//         #endif
//         #if DEBUG_PLOT_FF
//             I_ww[q] += w_scaling*w_scaling*water_water_sum*ff_table.index(neutral_oxygen_index, neutral_oxygen_index).evaluate(q);
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


ScatteringProfile CompositeDistanceHistogramFF::debye_transform() const {
    const auto& ff_table = form_factor::storage::get_precalculated_form_factor_table();
    const auto& sinqd_table = table::ArrayDebyeTable::get_default_table();

    // calculate the Debye scattering intensity
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    std::vector<double> Iq(debye_axis.bins, 0);
    std::vector<double> q_axis = debye_axis.as_vector();

    unsigned int ff_exv_index = static_cast<int>(form_factor::form_factor_t::EXCLUDED_VOLUME);
    unsigned int ff_water_index = static_cast<int>(form_factor::form_factor_t::O);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            // atom-atom
            for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                double atom_atom_sum = std::inner_product(cp_pp.begin(ff1, ff2), cp_pp.end(ff1, ff2), sinqd_table.begin(q), 0.0);
                Iq[q] += atom_atom_sum*ff_table.index(ff1, ff2).evaluate(q);
            }

            // atom-exv
            double atom_exv_sum = std::inner_product(cp_pp.begin(ff1, ff_exv_index), cp_pp.end(ff1, ff_exv_index), sinqd_table.begin(q), 0.0);
            Iq[q] -= 2*exv_scaling*atom_exv_sum*ff_table.index(ff1, ff_exv_index).evaluate(q);

            // atom-water
            double atom_water_sum = std::inner_product(cp_hp.begin(ff1), cp_hp.end(ff1), sinqd_table.begin(q), 0.0);
            Iq[q] += 2*w_scaling*atom_water_sum*ff_table.index(ff1, ff_water_index).evaluate(q);
        }

        // exv-exv
        double exv_exv_sum = std::inner_product(cp_pp.begin(ff_exv_index, ff_exv_index), cp_pp.end(ff_exv_index, ff_exv_index), sinqd_table.begin(q), 0.0);
        Iq[q] += exv_scaling*exv_scaling*exv_exv_sum*ff_table.index(ff_exv_index, ff_exv_index).evaluate(q);

        // exv-water
        double exv_water_sum = std::inner_product(cp_hp.begin(ff_exv_index), cp_hp.end(ff_exv_index), sinqd_table.begin(q), 0.0);
        Iq[q] -= 2*exv_scaling*w_scaling*exv_water_sum*ff_table.index(ff_exv_index, ff_water_index).evaluate(q);

        // water-water
        double water_water_sum = std::inner_product(cp_hh.begin(), cp_hh.end(), sinqd_table.begin(q), 0.0);
        Iq[q] += w_scaling*w_scaling*water_water_sum*ff_table.index(ff_water_index, ff_water_index).evaluate(q);
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
            std::transform(p_pp.begin(), p_pp.end(), cp_pp.begin(ff1, ff2), p_pp.begin(), std::plus<double>());
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