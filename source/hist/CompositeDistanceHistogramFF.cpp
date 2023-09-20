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
    : CompositeDistanceHistogram(std::vector<double>(0), std::vector<double>(0), std::vector<double>(0), std::move(p_tot), axis), p_pp(std::move(p_pp)), p_hp(std::move(p_hp)), p_hh(std::move(p_hh)) {}

CompositeDistanceHistogramFF::~CompositeDistanceHistogramFF() = default;

ScatteringHistogram CompositeDistanceHistogramFF::debye_transform() const {
    static Container2D<hist::detail::PrecalculatedFormFactorProduct> ff_table = hist::detail::PrecalculatedFormFactorProduct::generate_table();

    // calculate the Debye scattering intensity
    Axis debye_axis(settings::axes::qmin, settings::axes::qmax, settings::axes::bins);

    std::vector<double> Iq(debye_axis.bins, 0);
    std::vector<double> q_axis = debye_axis.as_vector();

    unsigned int excluded_volume_index = static_cast<int>(hist::detail::form_factor_t::EXCLUDED_VOLUME);
    unsigned int neutral_oxygen_index = static_cast<int>(hist::detail::form_factor_t::NEUTRAL_OXYGEN);
    for (unsigned int q = 0; q < debye_axis.bins; ++q) {
        for (unsigned int ff1 = 0; ff1 < detail::FormFactor::get_count()-1; ++ff1) {
            // atom-atom
            for (unsigned int ff2 = 0; ff2 < detail::FormFactor::get_count()-1; ++ff2) {
                double atom_atom_sum = 0;
                for (unsigned int d = 0; d < axis.bins; ++d) {
                    atom_atom_sum += p_pp.index(ff1, ff2, d)*sinqd_table->lookup(q, d);
                }
                Iq[q] += atom_atom_sum*ff_table.index(ff1, ff2).evaluate(q);
            }

            double atom_exv_sum = 0; // atom-exv
            double exv_exv_sum = 0;  // exv-exv
            for (unsigned int d = 0; d < axis.bins; ++d) {
                double sinqd = sinqd_table->lookup(q, d);
                atom_exv_sum += (p_pp.index(ff1, excluded_volume_index, d) + p_pp.index(excluded_volume_index, ff1, d))*sinqd;
                exv_exv_sum += p_pp.index(excluded_volume_index, excluded_volume_index, d)*sinqd;
            }
            Iq[q] -= atom_exv_sum*ff_table.index(ff1, excluded_volume_index).evaluate(q);
            Iq[q] += exv_exv_sum*ff_table.index(excluded_volume_index, excluded_volume_index).evaluate(q);

            // atom-water
            double val = 0;
            for (unsigned int d = 0; d < axis.bins; ++d) {
                val += p_hp.index(ff1, d)*sinqd_table->lookup(q, d);
            }
            Iq[q] += k*val*ff_table.index(ff1, neutral_oxygen_index).evaluate(q);
        }

        // exv-water
        double exv_water_sum = 0;
        for (unsigned int d = 0; d < axis.bins; ++d) {
            exv_water_sum += p_hp.index(excluded_volume_index, d)*sinqd_table->lookup(q, d);
        }
        Iq[q] -= k*2*exv_water_sum*ff_table.index(excluded_volume_index, neutral_oxygen_index).evaluate(q);

        // water-water
        double water_water_sum = 0;
        for (unsigned int d = 0; d < axis.bins; ++d) {
            water_water_sum += p_hh.index(d)*sinqd_table->lookup(q, d);
        }
        Iq[q] += k*k*water_water_sum*ff_table.index(neutral_oxygen_index, neutral_oxygen_index).evaluate(q);
    }

    // // atom-atom
    // for (unsigned int q = 0; q < debye_axis.bins; ++q) {
    //     for (unsigned int ff1 = 0; ff1 < detail::FormFactor::get_count(); ++ff1) {
    //         for (unsigned int ff2 = 0; ff2 < detail::FormFactor::get_count(); ++ff2) {
    //             for (unsigned int d = 0; d < axis.bins; ++d) {
    //                 Iq[q] += (
    //                     p_pp.index(ff1, ff2, d)*ff_table.index(ff1, ff2).evaluate(q) // Z1*Z2*f1(q)*f2(q)
    //                     + 5.95*5.95*std::exp(-1.62*1.62*q*q))                        // Zexv^2*ffexv^2
                    
    //                 *sinqd_table->lookup(q, d); 
    //             }
    //         }
    //     }
    // }

    // // atom-water
    // for (unsigned int q = 0; q < debye_axis.bins; ++q) {
    //     for (unsigned int ff1 = 0; ff1 < detail::FormFactor::get_count(); ++ff1) {
    //         for (unsigned int d = 0; d < axis.bins; ++d) {
    //             Iq[q] += (k*p_hp.index(ff1, d)*ff_table.index(ff1, static_cast<int>(hist::detail::form_factor_t::NEUTRAL_OXYGEN)).evaluate(q) - 5.95*std::exp(-1.62*q*q/2))*sinqd_table->lookup(q, d);
    //         }
    //     }
    // }

    // // water-water
    // for (unsigned int q = 0; q < debye_axis.bins; ++q) {
    //     for (unsigned int d = 0; d < axis.bins; ++d) {
    //         Iq[q] += (k2*p_hh.index(d)*ff_table.index(static_cast<int>(hist::detail::form_factor_t::NEUTRAL_OXYGEN), static_cast<int>(hist::detail::form_factor_t::NEUTRAL_OXYGEN)).evaluate(q)  - 5.95*std::exp(-1.62*q*q/2))*sinqd_table->lookup(q, d);
    //     }
    // }

    return ScatteringHistogram(Iq, debye_axis);
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