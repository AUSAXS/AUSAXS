#include <hist/DistanceHistogram.h>
#include <hist/CompositeDistanceHistogram.h>
#include <hist/Histogram.h>
#include <table/DebyeTable.h>
#include <dataset/SimpleDataset.h>
#include <settings/HistogramSettings.h>

using namespace hist;

DistanceHistogram::DistanceHistogram() = default;

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

    table::DebyeTable::check_default(q_axis, d_axis);
}

ScatteringProfile DistanceHistogram::debye_transform() const {
    // calculate the Debye scattering intensity
    Axis debye_axis(settings::axes::qmin, settings::axes::qmax, settings::axes::bins);
    const auto& sinqd_table = table::DebyeTable::get_default_table();

    // calculate the scattering intensity based on the Debye equation
    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int q = 0; q < debye_axis.bins; ++q) { // iterate through all q values
        Iq[q] = std::inner_product(p.begin(), p.end(), sinqd_table.begin(q), 0.0);
        Iq[q] *= std::exp(-q_axis[q]*q_axis[q]); // form factor
    }
    return ScatteringProfile(Iq, debye_axis);
}

// SimpleDataset DistanceHistogram::debye_transform(const std::vector<double>& q) const {
//     // calculate the scattering intensity based on the Debye equation
//     std::vector<double> Iq(q.size(), 0);
//     for (unsigned int i = 0; i < q.size(); ++i) { // iterate through all q values
//         for (unsigned int j = 0; j < p.size(); ++j) { // iterate through the distance histogram
//             double qd = q[i]*d_axis[j];
//             if (qd < 1e-6) {Iq[i] += 1;}
//             else {Iq[i] += p[j]*std::sin(qd)/qd;}
//         }
//         Iq[i] *= std::exp(-q[i]*q[i]); // form factor
//     }
//     return SimpleDataset(q, Iq);
// }

const std::vector<double>& DistanceHistogram::get_d_axis() const {return d_axis;}

const std::vector<double>& DistanceHistogram::get_q_axis() const {return q_axis;}

const std::vector<double>& DistanceHistogram::get_total_counts() const {return get_counts();}

std::vector<double>& DistanceHistogram::get_total_counts() {return get_counts();}