#include <hist/ScatteringHistogram.h>
#include <utility/Settings.h>
#include <utility/Constants.h>
#include <utility/Axis.h>

#include <utility>
#include <math.h>
#include <algorithm>
#include <iostream>

using namespace hist;

void ScatteringHistogram::setup() {
    // calculate what distance each bin represents
    d = std::vector<double>(axis.bins, 0);
    double d_width = axis.width();

    // we use the middle of each bin as the d-value, except for the very first one which we fix as 0 since it primarily contains self-correlation terms
    for (unsigned int i = 1; i < axis.bins; i++) {
        d[i] = axis.min + d_width*(i+0.5);
    }    

    // prepare the q values for the intensity calculations
    Axis debye_axis = Axis(setting::axes::bins, setting::axes::qmin, setting::axes::qmax);
    q = std::vector<double>(debye_axis.bins);
    double debye_width = debye_axis.width();
    for (unsigned int i = 0; i < debye_axis.bins; i++) {
        q[i] = debye_axis.min + i*debye_width;
    }

    sinqd_table.initialize(q, d);
}

void ScatteringHistogram::apply_water_scaling_factor(const double& k) {
    double k2 = pow(k, 2);
    for (unsigned int i = 0; i < axis.bins; i++) {p[i] = p_pp[i] + k*p_hp[i] + k2*p_hh[i];} // p = p_tot, inherited from Histogram
}

Histogram ScatteringHistogram::plot_debye_scattering() const {
    SimpleDataset data = calc_debye_scattering_intensity();
    return Histogram(data.y(), Axis(q.size(), q[0], q.back()));
}

SimpleDataset ScatteringHistogram::calc_debye_scattering_intensity(const std::vector<double>& q) const {
    // calculate the scattering intensity based on the Debye equation
    std::vector<double> Iq(q.size(), 0);
    for (unsigned int i = 0; i < q.size(); i++) { // iterate through all q values
        for (unsigned int j = 0; j < p.size(); j++) { // iterate through the distance histogram
            double qd = q[i]*d[j];
            if (qd < 1e-6) {Iq[i] += 1;}
            else {Iq[i] += p[j]*std::sin(qd)/qd;}
        }
        Iq[i] *= exp(-q[i]*q[i]); // form factor
    }
    return SimpleDataset(q, Iq, "q", "I");
}

SimpleDataset ScatteringHistogram::calc_debye_scattering_intensity() const {
    // calculate the Debye scattering intensity
    Axis debye_axis = Axis(setting::axes::bins, setting::axes::qmin, setting::axes::qmax);

    // calculate the scattering intensity based on the Debye equation
    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int i = 0; i < debye_axis.bins; i++) { // iterate through all q values
        for (unsigned int j = 0; j < p.size(); j++) { // iterate through the distance histogram
            Iq[i] += p[j]*sinqd_table.lookup(i, j);
        }
        Iq[i] *= exp(-q[i]*q[i]); // form factor
    }
    return SimpleDataset(q, Iq, "q", "I");
}

Histogram ScatteringHistogram::plot_guinier_approx() const {
    std::vector<double> Iq = calc_guinier_approx().col("logI");
    return Histogram(Iq, Axis(setting::axes::bins, setting::axes::qmin, setting::axes::qmax));
}

double ScatteringHistogram::calc_guinier_gyration_ratio_squared() const {
    double num = 0, denom = 0;
    for (unsigned int i = 0; i < p.size(); i++) {
        num += p[i]*pow(d[i], 2);
        denom += 2*p[i];
    }
    return num/denom;
}

SimpleDataset ScatteringHistogram::calc_guinier_approx() const {
    double Rg2 = calc_guinier_gyration_ratio_squared();

    Axis debye_axis = Axis(setting::axes::bins, setting::axes::qmin, setting::axes::qmax);
    std::vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int i = 0; i < debye_axis.bins; i++) { // iterate through all q values
        Iq[i] = std::exp(-pow(q[i], 2)*Rg2/3);
    }

    return SimpleDataset(q, Iq, "q", "I");
}

ScatteringHistogram& ScatteringHistogram::operator=(const ScatteringHistogram& h) = default;
ScatteringHistogram& ScatteringHistogram::operator=(ScatteringHistogram&& h) = default;