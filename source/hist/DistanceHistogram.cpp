#include <hist/DistanceHistogram.h>
#include <hist/Histogram.h>

using namespace hist;

DistanceHistogram::DistanceHistogram(std::vector<double>&& p_tot, const Axis& axis) : p_tot(std::move(p_tot)), axis(axis) {}

DistanceHistogram::DistanceHistogram(CompositeDistanceHistogram&& cdh) : p_tot(std::move(cdh.get_total_histogram())), axis(cdh.get_axis()) {}

Histogram DistanceHistogram::debye_transform() const {
    // calculate the scattering intensity based on the Debye equation
    std::vector<double> Iq(q.size(), 0);
    for (unsigned int i = 0; i < q.size(); ++i) { // iterate through all q values
        for (unsigned int j = 0; j < p.size(); ++j) { // iterate through the distance histogram
            double qd = q[i]*d[j];
            if (qd < 1e-6) {Iq[i] += 1;}
            else {Iq[i] += p[j]*std::sin(qd)/qd;}
        }
        Iq[i] *= std::exp(-q[i]*q[i]); // form factor
    }
    return Histogram(q, Iq);
}

Histogram DistanceHistogram::debye_transform(const std::vector<double>& q) const {
    return Histogram();
}

std::vector<double>& DistanceHistogram::get_total_histogram() {return p_tot;}

std::vector<double> DistanceHistogram::get_axis_vector() const {return axis.as_vector();}

Axis& DistanceHistogram::get_axis() {return axis;}

CompositeDistanceHistogram::CompositeDistanceHistogram(std::vector<double>&& p_pp, std::vector<double>&& p_hh, std::vector<double>&& p_hp, std::vector<double>&& p_tot, const Axis& axis) 
    : DistanceHistogram(std::move(p_tot), axis), p_pp(std::move(p_pp)), p_hh(std::move(p_hh)), p_hp(std::move(p_hp)) {}

Histogram CompositeDistanceHistogram::debye_transform() const {
    return Histogram();
}

Histogram CompositeDistanceHistogram::debye_transform(const std::vector<double>& q) const {
    return Histogram();
}

std::vector<double>& CompositeDistanceHistogram::get_pp_histogram() {return p_pp;}

std::vector<double>& CompositeDistanceHistogram::get_hh_histogram() {return p_hh;}

std::vector<double>& CompositeDistanceHistogram::get_hp_histogram() {return p_hp;}

void CompositeDistanceHistogram::apply_water_scaling_factor(double k) {
    double k2 = std::pow(k, 2);
    std::vector<double>& p_tot = get_total_histogram();
    for (unsigned int i = 0; i < get_axis().bins; ++i) {p_tot[i] = p_pp[i] + k*p_hp[i] + k2*p_hh[i];} // p = p_tot, inherited from Histogram
}

Histogram CompositeDistanceHistogramFF::debye_transform() const {
    return Histogram();
}

Histogram CompositeDistanceHistogramFF::debye_transform(const std::vector<double>& q) const {
    return Histogram();
}