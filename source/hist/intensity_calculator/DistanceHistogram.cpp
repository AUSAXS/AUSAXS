#include <hist/intensity_calculator/DistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <hist/Histogram.h>
#include <table/ArrayDebyeTable.h>
#include <dataset/SimpleDataset.h>
#include <settings/HistogramSettings.h>
#include <constants/Constants.h>

using namespace hist;

DistanceHistogram::DistanceHistogram() = default;

DistanceHistogram::DistanceHistogram(std::vector<constants::axes::d_type>&& p_tot, const Axis& axis) : Histogram(std::move(p_tot), axis) {
    initialize();
}

DistanceHistogram::DistanceHistogram(hist::Distribution1D&& p_tot, const Axis& axis) : DistanceHistogram(std::move(p_tot.get_container().get_data()), axis) {}

DistanceHistogram::DistanceHistogram(std::unique_ptr<ICompositeDistanceHistogram> cdh) : Histogram(std::move(cdh->get_total_counts()), cdh->get_axis()) {
    initialize();
}

DistanceHistogram::~DistanceHistogram() = default;

void DistanceHistogram::initialize() {
    q_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax).as_vector();
    d_axis = axis.as_vector();
    d_axis[0] = 0; // fix the first bin to 0 since it primarily contains self-correlation terms
    table::ArrayDebyeTable::check_default(q_axis, d_axis);
}

ScatteringProfile DistanceHistogram::debye_transform() const {
    // calculate the Debye scattering intensity
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    const auto& sinqd_table = table::ArrayDebyeTable::get_default_table();

    // calculate the scattering intensity based on the Debye equation
    std::vector<double> Iq(debye_axis.bins, 0);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) { // iterate through all q values
        Iq[q] = std::inner_product(p.begin(), p.end(), sinqd_table.begin(q), 0.0);
        Iq[q] *= std::exp(-q_axis[q]*q_axis[q]); // form factor
    }
    return ScatteringProfile(Iq, debye_axis);
}

const std::vector<double>& DistanceHistogram::get_d_axis() const {return d_axis;}

const std::vector<double>& DistanceHistogram::get_q_axis() const {return q_axis;}

const std::vector<double>& DistanceHistogram::get_total_counts() const {return get_counts();}

std::vector<double>& DistanceHistogram::get_total_counts() {return get_counts();}