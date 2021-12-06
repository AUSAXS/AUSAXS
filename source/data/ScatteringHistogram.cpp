#include "ScatteringHistogram.h"
#include "settings.h"
#include "constants.h"

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
    d = vector<double>(axes[0], 0);
    double d_width = (double) (axes[2]-axes[1])/axes[0];
    for (int i = 0; i < axes[0]; i++) {
        d[i] = axes[1] + d_width*i;
    }    

    // prepare the q values for the intensity calculations
    const vector<double>& debye_axes = setting::axes::scattering_intensity_plot_axes;
    q = vector<double>(debye_axes[0]);
    double debye_width = (double) (debye_axes[2]-debye_axes[1])/debye_axes[0];
    for (int i = 0; i < debye_axes[0]; i++) {
        q[i] = debye_axes[1] + i*debye_width;
    }
}

void ScatteringHistogram::apply_water_scaling_factor(const double& k) {
    double k2 = pow(k, 2);
    for (int i = 0; i < axes[0]; i++) {p_tot[i] = p_pp[i] + k*p_hp[i] + k2*p_hh[i];}
}

vector<shared_ptr<TH1D>> ScatteringHistogram::plot_distance() const {
    vector<shared_ptr<TH1D>> hists = {
        std::make_shared<TH1D>("h_pp", "hist", axes[0], axes[1], axes[2]), 
        std::make_shared<TH1D>("h_hh", "hist", axes[0], axes[1], axes[2]), 
        std::make_shared<TH1D>("h_hp", "hist", axes[0], axes[1], axes[2]), 
        std::make_shared<TH1D>("h_tot", "hist", axes[0], axes[1], axes[2])
    };

    for (int i = 1; i < axes[0]; i++) {
        hists[0]->SetBinContent(i, p_pp[i-1]);
        hists[1]->SetBinContent(i, p_hh[i-1]);
        hists[2]->SetBinContent(i, p_hp[i-1]);
        hists[3]->SetBinContent(i, p_tot[i-1]);
    }

    return hists;
}

unique_ptr<TH1D> ScatteringHistogram::plot_debye_scattering() const {
    vector<double> Iq = calc_debye_scattering_intensity();
    const vector<double>& debye_axes = setting::axes::scattering_intensity_plot_axes;
    unique_ptr<TH1D> h = std::make_unique<TH1D>("hI_debye", "hist", debye_axes[0], debye_axes[1], debye_axes[2]);

    for (size_t i = 0; i < Iq.size(); i++) {
        // in ROOT histograms, bin 0 is an underflow bin, and n+1 is an overflow bin
        h->SetBinContent(i+1, Iq[i]);
    }
    return h;
}

vector<double> ScatteringHistogram::calc_debye_scattering_intensity() const {
    // calculate the Debye scattering intensity
    const vector<double>& debye_axes = setting::axes::scattering_intensity_plot_axes;

    // calculate the scattering intensity based on the Debye equation
    vector<double> Iq(debye_axes[0], 0);
    for (int i = 0; i < debye_axes[0]; i++) { // iterate through all q values
        for (size_t j = 0; j < p_tot.size(); j++) { // iterate through the distance histogram
            if (q[i]*d[j] < 1e-9) { // if qd is very close to zero, we fix sin(qd)/qd to 1
                Iq[i] += p_tot[j];
            } else {
                Iq[i] += p_tot[j]*sin(q[i]*d[j])/(q[i]*d[j]);
            }
        }
    }
    return Iq;
}

unique_ptr<TH1D> ScatteringHistogram::plot_guinier_approx() const {
    vector<double> Iq = calc_guinier_approx();

    const vector<double>& debye_axes = setting::axes::scattering_intensity_plot_axes;
    unique_ptr<TH1D> h = std::make_unique<TH1D>("hI_guinier", "hist", debye_axes[0], debye_axes[1], debye_axes[2]);

    for (int i = 0; i < debye_axes[0]; i++) {
        // in ROOT histograms, bin 0 is an underflow bin, and n+1 is an overflow bin
        h->SetBinContent(i+1, Iq[i]);
    }
    return h;
}

double ScatteringHistogram::calc_guinier_gyration_ratio_squared() const {
    double num = 0, denom = 0;
    for (size_t i = 0; i < p_tot.size(); i++) {
        num += p_tot[i]*pow(d[i], 2);
        denom += 2*p_tot[i];
    }
    return num/denom;
}

vector<double> ScatteringHistogram::calc_guinier_approx() const {
    double Rg2 = calc_guinier_gyration_ratio_squared();

    const vector<double>& debye_axes = setting::axes::scattering_intensity_plot_axes;
    vector<double> Iq(debye_axes[0], 0);
    double log_conv = log10(exp(1)); // we have to convert natural log to log10
    for (int i = 0; i < debye_axes[0]; i++) { // iterate through all q values
        Iq[i] = -pow(q[i], 2)*Rg2/3*log_conv;
    }

    return Iq;
}