#include "ScatteringHistogram.h"
#include "settings.h"
#include "constants.h"
#include "data/Axis.h"

#include <utility>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <stdexcept>

#include "TH1D.h"
#include "TCanvas.h"

using std::cout, std::endl;
using namespace ROOT;

void ScatteringHistogram::setup() {
    // calculate what distance each bin represents
    _d = vector<double>(axis.bins);
    double d_width = axis.width();
    for (int i = 0; i < axis.bins; i++) {
        _d[i] = axis.min + d_width*i;
    }    

    // prepare the q values for the intensity calculations
    const Axis& debye_axis = setting::axes::scattering_intensity_plot_axis;
    _q = vector<double>(debye_axis.bins);
    double debye_width = debye_axis.width();
    for (int i = 0; i < debye_axis.bins; i++) {
        _q[i] = debye_axis.min + i*debye_width;
    }
}

void ScatteringHistogram::apply_water_scaling_factor(const double& k) {
    double k2 = pow(k, 2);
    for (int i = 0; i < axis.bins; i++) {p[i] = p_pp[i] + k*p_hp[i] + k2*p_hh[i];} // p = p_tot, inherited from Histogram
}

vector<shared_ptr<TH1D>> ScatteringHistogram::plot_distance() const {
    vector<shared_ptr<TH1D>> hists = {
        std::make_shared<TH1D>("h_pp", "hist", axis.bins, axis.min, axis.max), 
        std::make_shared<TH1D>("h_hh", "hist", axis.bins, axis.min, axis.max), 
        std::make_shared<TH1D>("h_hp", "hist", axis.bins, axis.min, axis.max), 
        std::make_shared<TH1D>("h_tot", "hist", axis.bins, axis.min, axis.max)
    };

    for (int i = 1; i < axis.bins; i++) { 
        hists[0]->SetBinContent(i, p_pp[i-1]);
        hists[1]->SetBinContent(i, p_hh[i-1]);
        hists[2]->SetBinContent(i, p_hp[i-1]);
        hists[3]->SetBinContent(i, p_tot[i-1]);
    }

    return hists;
}

unique_ptr<TH1D> ScatteringHistogram::plot_debye_scattering() const {
    vector<double> Iq = calc_debye_scattering_intensity();
    const Axis& debye_axis = setting::axes::scattering_intensity_plot_axis;
    unique_ptr<TH1D> h = std::make_unique<TH1D>("hI_debye", "hist", debye_axis.bins, debye_axis.min, debye_axis.max);

    for (unsigned int i = 0; i < Iq.size(); i++) {
        // in ROOT histograms, bin 0 is an underflow bin, and n+1 is an overflow bin
        std::cout << "Bin " << i << ": " << Iq[i] << std::endl;
        h->SetBinContent(i+1, Iq[i]);
    }
    return h;
}

vector<double> ScatteringHistogram::calc_debye_scattering_intensity() const {
    // calculate the Debye scattering intensity
    const Axis& debye_axis = setting::axes::scattering_intensity_plot_axis;

    // calculate the scattering intensity based on the Debye equation
    vector<double> Iq(debye_axis.bins, 0);
    for (int i = 0; i < debye_axis.bins; i++) { // iterate through all q values
        for (unsigned int j = 0; j < p_tot.size(); j++) { // iterate through the distance histogram
            if (q[i]*_d[j] < 1e-9) { // if qd is very close to zero, we fix sin(qd)/qd to 1
                Iq[i] += p_tot[j];
            } else {
                Iq[i] += p_tot[j]*sin(q[i]*_d[j])/(q[i]*_d[j]);
            }
        }
        Iq[i] *= exp(-q[i]*q[i]); // form factor
    }
    return Iq;
}

unique_ptr<TH1D> ScatteringHistogram::plot_guinier_approx() const {
    vector<double> Iq = calc_guinier_approx();

    const Axis& debye_axis = setting::axes::scattering_intensity_plot_axis;
    unique_ptr<TH1D> h = std::make_unique<TH1D>("hI_guinier", "hist", debye_axis.bins, debye_axis.min, debye_axis.max);

    for (int i = 0; i < debye_axis.bins; i++) {
        // in ROOT histograms, bin 0 is an underflow bin, and n+1 is an overflow bin
        h->SetBinContent(i+1, Iq[i]);
    }
    return h;
}

double ScatteringHistogram::calc_guinier_gyration_ratio_squared() const {
    double num = 0, denom = 0;
    for (unsigned int i = 0; i < p_tot.size(); i++) {
        num += p_tot[i]*pow(_d[i], 2);
        denom += 2*p_tot[i];
    }
    return num/denom;
}

vector<double> ScatteringHistogram::calc_guinier_approx() const {
    double Rg2 = calc_guinier_gyration_ratio_squared();

    const Axis& debye_axis = setting::axes::scattering_intensity_plot_axis;
    vector<double> Iq(debye_axis.bins, 0);
    double log_conv = log10(exp(1)); // we have to convert natural log to log10
    for (int i = 0; i < debye_axis.bins; i++) { // iterate through all q values
        Iq[i] = -pow(q[i], 2)*Rg2/3*log_conv;
    }

    return Iq;
}

ScatteringHistogram& ScatteringHistogram::operator=(const ScatteringHistogram& h) {
    return operator=(ScatteringHistogram(h));
}

ScatteringHistogram& ScatteringHistogram::operator=(ScatteringHistogram&& h) {
    p = std::move(h.p);
    _p_pp = std::move(h.p_pp);
    _p_hh = std::move(h.p_hh);
    _p_hp = std::move(h.p_hp);
    _q = std::move(h.q);
    _d = std::move(h._d);
    return *this;
}