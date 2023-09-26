#include <hist/CompositeDistanceHistogramFF.h>
#include <hist/CompositeDistanceHistogram.h>
#include <hist/Histogram.h>
#include <hist/DebyeLookupTable.h>
#include <hist/detail/FormFactor.h>
#include <hist/detail/PrecalculatedFormFactorProduct.h>
#include <dataset/SimpleDataset.h>
#include <data/Protein.h>
#include <data/Atom.h>
#include <settings/HistogramSettings.h>

using namespace hist;

CompositeDistanceHistogramFF::CompositeDistanceHistogramFF(Container3D<double>&& p_pp, Container2D<double>&& p_hp, Container1D<double>&& p_hh, std::vector<double>&& p_tot, const Axis& axis) 
    : CompositeDistanceHistogram(std::move(p_tot), axis), cp_pp(std::move(p_pp)), cp_hp(std::move(p_hp)), cp_hh(std::move(p_hh)) {}

CompositeDistanceHistogramFF::CompositeDistanceHistogramFF(CompositeDistanceHistogramFF&& other) noexcept 
    : CompositeDistanceHistogram(std::move(other.p_pp), std::move(other.p_hp), std::move(other.p_hh), std::move(other.p), other.axis), cp_pp(std::move(other.cp_pp)), cp_hp(std::move(other.cp_hp)), cp_hh(std::move(other.cp_hh)), k(other.k) {}

CompositeDistanceHistogramFF::~CompositeDistanceHistogramFF() = default;

#include <plots/PlotDataset.h>
#include <settings/GeneralSettings.h>
#define DEBUG_DEBYE_TRANSFORM 0
#if DEBUG_DEBYE_TRANSFORM
    unsigned int qcheck = 0;
#endif

ScatteringHistogram CompositeDistanceHistogramFF::debye_transform() const {
    static Container2D<hist::detail::PrecalculatedFormFactorProduct> ff_table = hist::detail::PrecalculatedFormFactorProduct::generate_table();

    // calculate the Debye scattering intensity
    Axis debye_axis(settings::axes::qmin, settings::axes::qmax, settings::axes::bins);

    std::vector<double> Iq(debye_axis.bins, 0);
    std::vector<double> q_axis = debye_axis.as_vector();

    std::vector<double> ff_effective(q_axis.size(), 0);

    unsigned int excluded_volume_index = static_cast<int>(hist::detail::form_factor_t::EXCLUDED_VOLUME);
    unsigned int neutral_oxygen_index = static_cast<int>(hist::detail::form_factor_t::NEUTRAL_OXYGEN);
    for (unsigned int q = 0; q < debye_axis.bins; ++q) {
        for (unsigned int ff1 = 0; ff1 < detail::FormFactor::get_count_without_excluded_volume(); ++ff1) {
            // atom-atom
            for (unsigned int ff2 = 0; ff2 < detail::FormFactor::get_count_without_excluded_volume(); ++ff2) {
                double atom_atom_sum = 0;
                for (unsigned int d = 0; d < axis.bins; ++d) {
                    atom_atom_sum += cp_pp.index(ff1, ff2, d)*sinqd_table->lookup(q, d);
                }
                Iq[q] += atom_atom_sum*ff_table.index(ff1, ff2).evaluate(q);

                #if DEBUG_DEBYE_TRANSFORM
                    if (q==qcheck && atom_atom_sum != 0) {
                        std::cout << "(aa) Iq[" << q << "] += " << atom_atom_sum << " * ff_table[" << ff1 << ", " << ff2 << "](" << q << ") = " << atom_atom_sum*ff_table.index(ff1, ff2).evaluate(q) << std::endl;
                        for (unsigned int d = 0; d < axis.bins; ++d) {
                            if (cp_pp.index(ff1, ff2, d) != 0) {std::cout << "\taa_sum += p_pp[" << ff1 << ", " << ff2 << ", " << d << "] = " << cp_pp.index(ff1, ff2, d) << "*" << sinqd_table->lookup(q, d) << std::endl;}
                        }
                        std::cout << "\taasum = " << atom_atom_sum << std::endl;
                    }
                #endif
                ff_effective[q] += atom_atom_sum*ff_table.index(ff1, ff2).evaluate(q);
            }

            double atom_exv_sum = 0; // atom-exv
            for (unsigned int d = 0; d < axis.bins; ++d) {
                //! factor 2 or not?
                atom_exv_sum += 2*(cp_pp.index(ff1, excluded_volume_index, d) + cp_pp.index(excluded_volume_index, ff1, d))*sinqd_table->lookup(q, d);
            }
            Iq[q] -= atom_exv_sum*ff_table.index(ff1, excluded_volume_index).evaluate(q);

            #if DEBUG_DEBYE_TRANSFORM
                if (q==qcheck && atom_exv_sum != 0) {
                    std::cout << "(ae) Iq[" << q << "] -= " << atom_exv_sum << " * ff_table[" << ff1 << ", " << excluded_volume_index << "](" << q << ") = " << atom_exv_sum*ff_table.index(ff1, excluded_volume_index).evaluate(q) << std::endl;
                    for (unsigned int d = 0; d < axis.bins; ++d) {
                        if (atom_exv_sum != 0) {std::cout << "\tae_sum += 2*p_pp[" << ff1 << ", " << excluded_volume_index << ", " << d << "] = " << 2*(cp_pp.index(ff1, excluded_volume_index, d) + cp_pp.index(excluded_volume_index, ff1, d)) << "*" << sinqd_table->lookup(q, d) << std::endl;}
                    }
                    std::cout << "\taesum = " << atom_exv_sum << std::endl;
                }
            #endif

            // atom-water
            double atom_water_sum = 0;
            for (unsigned int d = 0; d < axis.bins; ++d) {
                atom_water_sum += cp_hp.index(ff1, d)*sinqd_table->lookup(q, d);
            }
            Iq[q] += k*atom_water_sum*ff_table.index(ff1, neutral_oxygen_index).evaluate(q);

            #if DEBUG_DEBYE_TRANSFORM
                if (q==qcheck && ff1 == 1) {
                    std::cout << "(aw) Iq[" << q << "] += " << k*atom_water_sum << " * ff_table[" << ff1 << ", " << neutral_oxygen_index << "](" << q << ") = " << k*atom_water_sum*ff_table.index(ff1, neutral_oxygen_index).evaluate(q) << std::endl;
                    for (unsigned int d = 0; d < axis.bins; ++d) {
                        if (atom_water_sum != 0) {std::cout << "\taw_sum += p_hp[" << ff1 << ", " << d << "] = " << cp_hp.index(ff1, d) << "*" << sinqd_table->lookup(q, d) << " = " << cp_hp.index(ff1, d)*sinqd_table->lookup(q, d) << std::endl;}
                    }
                    std::cout << "\tawsum = " << atom_water_sum << std::endl;
                }
            #endif
        }

        // exv-exv
        double exv_exv_sum = 0;
        for (unsigned int d = 0; d < axis.bins; ++d) {
            exv_exv_sum += cp_pp.index(excluded_volume_index, excluded_volume_index, d)*sinqd_table->lookup(q, d);
        }
        Iq[q] += exv_exv_sum*ff_table.index(excluded_volume_index, excluded_volume_index).evaluate(q);
        ff_effective[q] /= exv_exv_sum*ff_table.index(excluded_volume_index, excluded_volume_index).evaluate(q);

        #if DEBUG_DEBYE_TRANSFORM
            if (q==qcheck) {
                std::cout << "(ee) Iq[" << q << "] += " << exv_exv_sum << " * ff_table[" << excluded_volume_index << ", " << excluded_volume_index << "](" << q << ") = " << exv_exv_sum*ff_table.index(excluded_volume_index, excluded_volume_index).evaluate(q) << std::endl;
                for (unsigned int d = 0; d < axis.bins; ++d) {
                    if (exv_exv_sum != 0) {std::cout << "\tee_sum += p_pp[" << excluded_volume_index << ", " << excluded_volume_index << ", " << d << "] = " << cp_pp.index(excluded_volume_index, excluded_volume_index, d) << std::endl;}
                }
                std::cout << "\teesum = " << exv_exv_sum << std::endl;
            }
        #endif

        // exv-water
        double exv_water_sum = 0;
        for (unsigned int d = 0; d < axis.bins; ++d) {
            exv_water_sum += cp_hp.index(excluded_volume_index, d)*sinqd_table->lookup(q, d);
        }
        Iq[q] -= k*2*exv_water_sum*ff_table.index(excluded_volume_index, neutral_oxygen_index).evaluate(q);

        #if DEBUG_DEBYE_TRANSFORM
            if (q==qcheck) {
                std::cout << "(ew) Iq[" << q << "] -= " << k*2*exv_water_sum << " * ff_table[" << excluded_volume_index << ", " << neutral_oxygen_index << "](" << q << ") = " << k*2*exv_water_sum*ff_table.index(excluded_volume_index, neutral_oxygen_index).evaluate(q) << std::endl;
                for (unsigned int d = 0; d < axis.bins; ++d) {
                    std::cout << "\tew_sum += p_hp[" << excluded_volume_index << ", " << d << "] = " << cp_hp.index(excluded_volume_index, d) << "*" << sinqd_table->lookup(q, d) << std::endl;
                }
                std::cout << "\tewsum = " << exv_water_sum << std::endl;
            }
        #endif

        // water-water
        double water_water_sum = 0;
        for (unsigned int d = 0; d < axis.bins; ++d) {
            water_water_sum += cp_hh.index(d)*sinqd_table->lookup(q, d);
        }
        Iq[q] += k*k*water_water_sum*ff_table.index(neutral_oxygen_index, neutral_oxygen_index).evaluate(q);

        #if DEBUG_DEBYE_TRANSFORM
            if (q==qcheck) {std::cout << "(ww) Iq[" << q << "] += " << k*k*water_water_sum << " * ff_table[" << neutral_oxygen_index << ", " << neutral_oxygen_index << "](" << q << ") = " << k*k*water_water_sum*ff_table.index(neutral_oxygen_index, neutral_oxygen_index).evaluate(q) << std::endl;}
        #endif
    }
    SimpleDataset temp(q_axis, ff_effective);
    temp.add_plot_options({{"xlabel", "q"}, {"ylabel", "ff_effective"}});
    plots::PlotDataset::quick_plot(temp, settings::general::output + "ff_effective.png");
    return ScatteringHistogram(Iq, debye_axis);
}

const std::vector<double>& CompositeDistanceHistogramFF::get_counts() const {
    p = std::vector<double>(axis.bins, 0);
    auto& p_pp = get_pp_counts();
    auto& p_hp = get_hp_counts();
    auto& p_hh = get_hh_counts();
    for (unsigned int i = 0; i < axis.bins; ++i) {
        p[i] = p_pp[i] + k*p_hp[i] + k*k*p_hh[i];
    }
    return p.data;
}

const std::vector<double>& CompositeDistanceHistogramFF::get_pp_counts() const {
    p_pp = std::vector<double>(axis.bins, 0);
    for (unsigned int i = 0; i < axis.bins; ++i) {
        for (unsigned int ff1 = 0; ff1 < hist::detail::FormFactor::get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < hist::detail::FormFactor::get_count_without_excluded_volume(); ++ff2) {
                p_pp[i] += cp_pp.index(ff1, ff2, i);
            }
        }
    }
    return p_pp;
}

const std::vector<double>& CompositeDistanceHistogramFF::get_hp_counts() const {
    p_hp = std::vector<double>(axis.bins, 0);
    for (unsigned int i = 0; i < axis.bins; ++i) {
        for (unsigned int ff1 = 0; ff1 < hist::detail::FormFactor::get_count_without_excluded_volume(); ++ff1) {
            p_hp[i] += cp_hp.index(ff1, i);
        }
    }
    return p_hp;
}

const std::vector<double>& CompositeDistanceHistogramFF::get_hh_counts() const {
    p_hh = std::vector<double>(axis.bins, 0);
    for (unsigned int i = 0; i < axis.bins; ++i) {
        p_hh[i] += cp_hh.index(i);
    }
    return p_hh;
}

SimpleDataset CompositeDistanceHistogramFF::debye_transform(const std::vector<double>& q) const {
    // calculate the scattering intensity based on the Debye equation
    double k2 = k*k;
    throw std::runtime_error("Not implemented");
    // std::vector<double> Iq(q.size(), 0);
    // for (unsigned int ff1 = 0; ff1 < detail::FormFactor::get_count(); ++ff1) {
    //     for (unsigned int ff2 = 0; ff2 < detail::FormFactor::get_count(); ++ff2) {
    //         auto ff1_func = detail::FormFactorStorage::get_form_factor(static_cast<detail::form_factor_t>(ff1));
    //         auto ff2_func = detail::FormFactorStorage::get_form_factor(static_cast<detail::form_factor_t>(ff2));
    //         for (unsigned int i = 0; i < q.size(); ++i) { // iterate through all q values
    //             for (unsigned int j = 0; j < axis.bins; ++j) { // iterate through the distance histogram
    //                 Iq[i] += p_pp.index(ff1, ff2, j)*sinqd_table->lookup(i, j);
    //             }
    //             Iq[i] *= ff1_func.evaluate(q[i])*ff2_func.evaluate(q[i]); // form factor
    //         }
    //     }
    // }
    // return SimpleDataset(Iq, q);
}

void CompositeDistanceHistogramFF::apply_water_scaling_factor(double k) {
    this->k = k;
}