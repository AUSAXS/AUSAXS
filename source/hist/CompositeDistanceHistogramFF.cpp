#include <hist/CompositeDistanceHistogramFF.h>
#include <hist/CompositeDistanceHistogram.h>
#include <hist/Histogram.h>
#include <hist/DebyeLookupTable.h>
#include <hist/detail/FormFactor.h>
#include <hist/detail/PrecalculatedFormFactorProduct.h>
#include <dataset/SimpleDataset.h>
#include <settings/HistogramSettings.h>

using namespace hist;

CompositeDistanceHistogramFF::CompositeDistanceHistogramFF(Container3D<double>&& p_pp, Container2D<double>&& p_hp, Container1D<double>&& p_hh, std::vector<double>&& p_tot, const Axis& axis) 
    : CompositeDistanceHistogram(std::vector<double>(0), std::vector<double>(0), std::vector<double>(0), std::move(p_tot), axis), p_pp(std::move(p_pp)), p_hp(std::move(p_hp)), p_hh(std::move(p_hh)) {}

CompositeDistanceHistogramFF::~CompositeDistanceHistogramFF() = default;

ScatteringHistogram CompositeDistanceHistogramFF::debye_transform() const {
    static Container2D<hist::detail::PrecalculatedFormFactorProduct> ff_table = hist::detail::PrecalculatedFormFactorProduct::generate_table();

    // calculate the Debye scattering intensity
    Axis debye_axis(settings::axes::qmin, settings::axes::qmax, settings::axes::bins);

    // calculate the scattering intensity based on the Debye equation
    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int ff1 = 0; ff1 < detail::FormFactor::get_count(); ++ff1) {
        for (unsigned int ff2 = 0; ff2 < detail::FormFactor::get_count(); ++ff2) {
            auto precalculated_ff = ff_table.index(ff1, ff2);
            for (unsigned int i = 0; i < debye_axis.bins; ++i) { // iterate through all q values
                for (unsigned int j = 0; j < axis.bins; ++j) { // iterate through the distance histogram
                    Iq[i] += p_pp.index(ff1, ff2, j)*sinqd_table->lookup(i, j);
                }
                Iq[i] *= precalculated_ff(i); // form factor
            }
        }
    }
    return ScatteringHistogram(Iq, debye_axis);
}

SimpleDataset CompositeDistanceHistogramFF::debye_transform(const std::vector<double>& q) const {
    // calculate the scattering intensity based on the Debye equation
    std::vector<double> Iq(q.size(), 0);
    for (unsigned int ff1 = 0; ff1 < detail::FormFactor::get_count(); ++ff1) {
        for (unsigned int ff2 = 0; ff2 < detail::FormFactor::get_count(); ++ff2) {
            auto ff1_func = detail::FormFactorStorage::get_form_factor(static_cast<detail::form_factor_t>(ff1));
            auto ff2_func = detail::FormFactorStorage::get_form_factor(static_cast<detail::form_factor_t>(ff2));
            for (unsigned int i = 0; i < q.size(); ++i) { // iterate through all q values
                for (unsigned int j = 0; j < axis.bins; ++j) { // iterate through the distance histogram
                    Iq[i] += p_pp.index(ff1, ff2, j)*sinqd_table->lookup(i, j);
                }
                Iq[i] *= ff1_func.evaluate(q[i])*ff2_func.evaluate(q[i]); // form factor
            }
        }
    }
    return SimpleDataset(Iq, q);
}

void CompositeDistanceHistogramFF::apply_water_scaling_factor(double k) {
    double k2 = std::pow(k, 2);
    auto& p_tot = get_total_histogram();
    for (unsigned int i = 0; i < axis.bins; ++i) {
        p_tot[i] = 0;
        for (unsigned int ff1 = 0; ff1 < detail::FormFactor::get_count(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < detail::FormFactor::get_count(); ++ff2) {
                p_tot[i] += p_pp.index(ff1, ff2, i);
            }
            p_tot[i] += k*p_hp.index(ff1, i);
        }
        p_tot[i] += k2*p_hh.index(i);
    }
}