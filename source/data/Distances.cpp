#include "Distances.h"
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

unique_ptr<TH1D> Distances::fit_debye_plot() const {
    // calculate the scattering intensity based on the Debye equation
    double I0;
    double r2 = pow(constants::radius::electron, 2);
    for (int j = 0; j < axes[0]; j++) { // iterate through the distance histogram
        I0 += binned_tot[j]*r2;
    }

    const double M = protein->get_mass();
    I0 /= M;

    vector<double> Iq = calc_debye_scattering_intensity();
    cout << "c is " << I0/(Iq[0]*r2);

    const vector<double>& debye_axes = setting::axes::scattering_intensity_plot_axes;
    unique_ptr<TH1D> h = std::make_unique<TH1D>("debye_fit", "hist", debye_axes[0], debye_axes[1], debye_axes[2]);
    return std::move(h);
}

vector<shared_ptr<TH1D>> Distances::plot_distance() const {
    vector<shared_ptr<TH1D>> hists = {
        std::make_shared<TH1D>("h_pp", "hist", axes[0], axes[1], axes[2]), 
        std::make_shared<TH1D>("h_hh", "hist", axes[0], axes[1], axes[2]), 
        std::make_shared<TH1D>("h_hp", "hist", axes[0], axes[1], axes[2]), 
        std::make_shared<TH1D>("h_tot", "hist", axes[0], axes[1], axes[2])
    };

    for (int i = 1; i < axes[0]; i++) {
        hists[0]->SetBinContent(i, binned_pp[i-1]);
        hists[1]->SetBinContent(i, binned_hh[i-1]);
        hists[2]->SetBinContent(i, binned_hp[i-1]);
        hists[3]->SetBinContent(i, binned_tot[i-1]);
    }

    return hists;
}

unique_ptr<TH1D> Distances::plot_debye_scattering() const {
    vector<double> Iq = calc_debye_scattering_intensity();
    const vector<double>& debye_axes = setting::axes::scattering_intensity_plot_axes;
    unique_ptr<TH1D> h = std::make_unique<TH1D>("hI_debye", "hist", debye_axes[0], debye_axes[1], debye_axes[2]);

    for (int i = 0; i < Iq.size(); i++) {
        // in ROOT histograms, bin 0 is an underflow bin, and n+1 is an overflow bin
        h->SetBinContent(i+1, Iq[i]);
    }
    return std::move(h);
}

vector<double> Distances::calc_debye_scattering_intensity() const {
    // calculate the Debye scattering intensity
    const vector<double>& debye_axes = setting::axes::scattering_intensity_plot_axes;

    // calculate what distance each bin represents
    vector<double> d(axes[0], 0);
    double p_width = (double) (axes[2]-axes[1])/axes[0];
    for (int i = 0; i < axes[0]; i++) {
        d[i] = axes[1] + p_width*i;
    }

    // calculate the scattering intensity based on the Debye equation
    vector<double> Iq(debye_axes[0], 0);
    double debye_width = (double) (debye_axes[2]-debye_axes[1])/debye_axes[0];
    for (int i = 0; i < debye_axes[0]; i++) { // iterate through all q values
        double q = debye_axes[1] + i*debye_width; // set the q value for this iteration
        for (int j = 0; j < axes[0]; j++) { // iterate through the distance histogram
            if (q*d[j] < 1e-9) { // if qd is very close to zero, we fix sin(qd)/qd to 1
                Iq[i] += binned_tot[j];
            } else {
                Iq[i] += binned_tot[j]*sin(q*d[j])/(q*d[j]);
            }
        }
    }
    return Iq;
}

unique_ptr<TH1D> Distances::plot_guinier_approx() const {
    vector<double> Iq = calc_guinier_approx();

    const vector<double>& debye_axes = setting::axes::scattering_intensity_plot_axes;
    unique_ptr<TH1D> h = std::make_unique<TH1D>("hI_guinier", "hist", debye_axes[0], debye_axes[1], debye_axes[2]);

    for (int i = 0; i < debye_axes[0]; i++) {
        // in ROOT histograms, bin 0 is an underflow bin, and n+1 is an overflow bin
        h->SetBinContent(i+1, Iq[i]);
    }
    return std::move(h);
}

double Distances::calc_guinier_gyration_ratio_squared() const {
    // calculate what distance each bin represents
    vector<double> d(axes[0], 0);
    double d_width = (double) (axes[2]-axes[1])/axes[0];
    for (int i = 0; i < axes[0]; i++) {
        d[i] = axes[1] + d_width*i;
    }
    
    double num = 0, denom = 0;
    for (int i = 0; i < axes[0]; i++) {
        num += binned_tot[i]*pow(d[i], 2);
        denom += 2*binned_tot[i];
    }
    return num/denom;
}

vector<double> Distances::calc_guinier_approx() const {
    double Rg2 = calc_guinier_gyration_ratio_squared();

    const vector<double>& debye_axes = setting::axes::scattering_intensity_plot_axes;
    vector<double> Iq(debye_axes[0], 0);
    double debye_width = (double) (debye_axes[2]-debye_axes[1])/debye_axes[0];
    double log_conv = log10(exp(1)); // we have to convert natural log to log10
    for (int i = 0; i < debye_axes[0]; i++) { // iterate through all q values
        double q = debye_axes[1] + i*debye_width; // set the q value for this iteration
        Iq[i] = -pow(q, 2)*Rg2/3*log_conv;
    }

    return Iq;
}

vector<double> Distances::get_xaxis() const {
    const vector<double>& debye_axes = setting::axes::scattering_intensity_plot_axes;
    vector<double> x(debye_axes[0], 0);
    double debye_width = (double) (debye_axes[2]-debye_axes[1])/debye_axes[0];
    for (int i = 0; i < debye_axes[0]; i++) {
        x[i] = debye_axes[1] + i*debye_width;
    }
    return x;
}

void Distances::set_axes(const vector<int> axes) {
    this->axes = axes;
    bin(axes);
}

void Distances::bin(const vector<int>& axes) {
    vector<double> p_pp(axes[0], 0);
    vector<double> p_hh(axes[0], 0);
    vector<double> p_hp(axes[0], 0);
    vector<double> p_tot(axes[0], 0);
    double width = (double) (axes[2]-axes[1])/axes[0]; // very important to cast this operation to a double - divison by two ints
    for (int i = 0; i < d_pp.size(); i++) {
        p_pp[std::round(d_pp[i]/width)] += w_pp[i];
    }
    for (int i = 0; i < d_hh.size(); i++) {
        p_hh[std::round(d_hh[i]/width)] += w_hh[i];
    }
    for (int i = 0; i < d_hp.size(); i++) {
        p_hp[std::round(d_hp[i]/width)] += w_hp[i];
    }
    for (int i = 0; i < axes[0]; i++) {
        p_tot[i] = p_pp[i] + p_hh[i] + p_hp[i];
    }

    this->binned_pp = p_pp;
    this->binned_hh = p_hh;
    this->binned_hp = p_hp;
    this->binned_tot = p_tot;
}