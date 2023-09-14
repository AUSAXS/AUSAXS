#include <hist/CompositeDistanceHistogramFF.h>
#include <hist/CompositeDistanceHistogram.h>
#include <hist/Histogram.h>
#include <hist/detail/FormFactorType.h>
#include <hist/DebyeLookupTable.h>
#include <settings/HistogramSettings.h>

using namespace hist;

CompositeDistanceHistogramFF::CompositeDistanceHistogramFF(Container3D<double>&& p_pp, Container2D<double>&& p_hp, Container1D<double>&& p_hh, std::vector<double>&& p_tot, const Axis& axis) 
    : CompositeDistanceHistogram(std::vector<double>(0), std::vector<double>(0), std::vector<double>(0), std::move(p_tot), axis), p_pp(std::move(p_pp)), p_hp(std::move(p_hp)), p_hh(std::move(p_hh)) {}

CompositeDistanceHistogramFF::~CompositeDistanceHistogramFF() = default;

Histogram CompositeDistanceHistogramFF::debye_transform() const {
    // calculate the Debye scattering intensity
    Axis debye_axis(settings::axes::qmin, settings::axes::qmax, settings::axes::bins);

    // calculate the scattering intensity based on the Debye equation
    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int ff1 = 0; ff1 < detail::ff::get_count(); ++ff1) {
        for (unsigned int ff2 = 0; ff2 < detail::ff::get_count(); ++ff2) {
            for (unsigned int i = 0; i < debye_axis.bins; ++i) { // iterate through all q values
                for (unsigned int j = 0; j < axis.bins; ++j) { // iterate through the distance histogram
                    auto pp = p_pp.index(ff1, ff2, j);
                    Iq[i] += p_pp.index(ff1, ff2, j)*sinqd_table->lookup(i, j)*detail::ff::get_form_factor(ff1, ff2, q_axis[i]);
                }
                Iq[i] *= std::exp(-q_axis[i]*q_axis[i]); // form factor
            }
        }
    }

    for (unsigned int i = 0; i < debye_axis.bins; ++i) { // iterate through all q values
        for (unsigned int j = 0; j < p.size(); ++j) { // iterate through the distance histogram
            Iq[i] += p[j]*sinqd_table->lookup(i, j);
        }
        Iq[i] *= std::exp(-q_axis[i]*q_axis[i]); // form factor
    }
    return Histogram(Iq, debye_axis);
}

Histogram CompositeDistanceHistogramFF::debye_transform(const std::vector<double>& q) const {
    return Histogram();
}