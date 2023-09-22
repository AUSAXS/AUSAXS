#include <hist/DistanceHistogram.h>
#include <hist/CompositeDistanceHistogram.h>
#include <hist/Histogram.h>
#include <hist/DebyeLookupTable.h>
#include <dataset/SimpleDataset.h>
#include <settings/HistogramSettings.h>

using namespace hist;

DistanceHistogram::DistanceHistogram(std::vector<double>&& p_tot, const Axis& axis) : Histogram(std::move(p_tot), axis) {
    initialize();
}

DistanceHistogram::DistanceHistogram(CompositeDistanceHistogram&& cdh) : Histogram(std::move(cdh.get_total_counts()), cdh.get_axis()) {
    initialize();
}

DistanceHistogram::~DistanceHistogram() = default;

void DistanceHistogram::initialize() {
    q_axis = Axis(settings::axes::qmin, settings::axes::qmax, settings::axes::bins).as_vector();
    d_axis = axis.as_vector(0.5);
    d_axis[0] = 0; // fix the first bin to 0 since it primarily contains self-correlation terms

    sinqd_table = std::make_unique<table::DebyeLookupTable>(q_axis, d_axis);
}

ScatteringHistogram DistanceHistogram::debye_transform() const {
    // calculate the Debye scattering intensity
    Axis debye_axis(settings::axes::qmin, settings::axes::qmax, settings::axes::bins);

    // calculate the scattering intensity based on the Debye equation
    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int i = 0; i < debye_axis.bins; ++i) { // iterate through all q values
        for (unsigned int j = 0; j < p.size(); ++j) { // iterate through the distance histogram
            Iq[i] += p[j]*sinqd_table->lookup(i, j);
        }
        Iq[i] *= std::exp(-q_axis[i]*q_axis[i]); // form factor
    }
    return ScatteringHistogram(Iq, debye_axis);
}

SimpleDataset DistanceHistogram::debye_transform(const std::vector<double>& q) const {
    // calculate the scattering intensity based on the Debye equation
    std::vector<double> Iq(q.size(), 0);
    for (unsigned int i = 0; i < q.size(); ++i) { // iterate through all q values
        for (unsigned int j = 0; j < p.size(); ++j) { // iterate through the distance histogram
            double qd = q[i]*d_axis[j];
            if (qd < 1e-6) {Iq[i] += 1;}
            else {Iq[i] += p[j]*std::sin(qd)/qd;}
        }
        Iq[i] *= std::exp(-q[i]*q[i]); // form factor
    }
    return SimpleDataset(q, Iq);
}

const std::vector<double>& DistanceHistogram::get_d_axis() const {return d_axis;}

const std::vector<double>& DistanceHistogram::get_q_axis() const {return q_axis;}

const std::vector<double>& DistanceHistogram::get_total_counts() const {return get_counts();}

std::vector<double>& DistanceHistogram::get_total_counts() {return get_counts();}