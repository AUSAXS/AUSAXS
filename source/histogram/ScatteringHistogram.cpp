#include <histogram/ScatteringHistogram.h>
#include <utility/Settings.h>
#include <utility/Constants.h>
#include <utility/Axis.h>

#include <utility>
#include <math.h>
#include <algorithm>
#include <iostream>

#include <TH1D.h>
#include <TCanvas.h>

using std::vector;
using namespace ROOT;
using namespace hist;

void ScatteringHistogram::setup() {
    // calculate what distance each bin represents
    d = vector<double>(axis.bins, 0);
    double d_width = axis.width();

    // we use the middle of each bin as the d-value, except for the very first one which we fix as 0 since it primarily contains self-correlation terms
    for (unsigned int i = 1; i < axis.bins; i++) {
        d[i] = axis.min + d_width*(i+0.5);
    }    

    // prepare the q values for the intensity calculations
    const Axis& debye_axis = setting::axes::scattering_intensity_plot_axis;
    q = vector<double>(debye_axis.bins);
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

vector<std::shared_ptr<TH1D>> ScatteringHistogram::plot_distance() const {
    vector<std::shared_ptr<TH1D>> hists = {
        std::make_shared<TH1D>("h_pp", "hist", axis.bins, axis.min, axis.max), 
        std::make_shared<TH1D>("h_hh", "hist", axis.bins, axis.min, axis.max), 
        std::make_shared<TH1D>("h_hp", "hist", axis.bins, axis.min, axis.max), 
        std::make_shared<TH1D>("h_tot", "hist", axis.bins, axis.min, axis.max)
    };

    for (unsigned int i = 1; i < axis.bins; i++) { 
        hists[0]->SetBinContent(i, p_pp[i-1]);
        hists[1]->SetBinContent(i, p_hh[i-1]);
        hists[2]->SetBinContent(i, p_hp[i-1]);
        hists[3]->SetBinContent(i, p[i-1]);
    }

    return hists;
}

std::unique_ptr<TH1D> ScatteringHistogram::plot_debye_scattering() const {
    SAXSDataset data = calc_debye_scattering_intensity();
    vector<double> I = data.get("I");
    vector<double> q = data.get("q");
    std::unique_ptr<TH1D> h = std::make_unique<TH1D>("hI_debye", "hist", q.size(), q[0], q[q.size()-1]);

    for (unsigned int i = 0; i < I.size(); i++) {
        // in ROOT histograms, bin 0 is an underflow bin, and n+1 is an overflow bin
        // std::cout << "Bin " << i << ": " << I[i] << std::endl;
        h->SetBinContent(i+1, I[i]);
    }
    return h;
}

SAXSDataset ScatteringHistogram::calc_debye_scattering_intensity(const vector<double>& q) const {
    // calculate the scattering intensity based on the Debye equation
    vector<double> Iq(q.size(), 0);
    for (unsigned int i = 0; i < q.size(); i++) { // iterate through all q values
        for (unsigned int j = 0; j < p.size(); j++) { // iterate through the distance histogram
            double qd = q[i]*d[j];
            if (qd < 1e-6) {Iq[i] += 1;}
            else {Iq[i] += p[j]*std::sin(qd)/qd;}
        }
        Iq[i] *= exp(-q[i]*q[i]); // form factor
    }
    return SAXSDataset(q, Iq);
}

SAXSDataset ScatteringHistogram::calc_debye_scattering_intensity() const {
    // calculate the Debye scattering intensity
    const Axis& debye_axis = setting::axes::scattering_intensity_plot_axis;

    // calculate the scattering intensity based on the Debye equation
    vector<double> Iq(debye_axis.bins, 0);
    for (unsigned int i = 0; i < debye_axis.bins; i++) { // iterate through all q values
        for (unsigned int j = 0; j < p.size(); j++) { // iterate through the distance histogram
            Iq[i] += p[j]*sinqd_table.lookup(i, j);
        }
        Iq[i] *= exp(-q[i]*q[i]); // form factor
    }
    return SAXSDataset(q, Iq);
}

std::unique_ptr<TH1D> ScatteringHistogram::plot_guinier_approx() const {
    vector<double> Iq = calc_guinier_approx().get("logI");

    const Axis& debye_axis = setting::axes::scattering_intensity_plot_axis;
    std::unique_ptr<TH1D> h = std::make_unique<TH1D>("hI_guinier", "hist", debye_axis.bins, debye_axis.min, debye_axis.max);

    for (unsigned int i = 0; i < debye_axis.bins; i++) {
        // in ROOT histograms, bin 0 is an underflow bin, and n+1 is an overflow bin
        h->SetBinContent(i+1, Iq[i]);
    }
    return h;
}

double ScatteringHistogram::calc_guinier_gyration_ratio_squared() const {
    double num = 0, denom = 0;
    for (unsigned int i = 0; i < p.size(); i++) {
        num += p[i]*pow(d[i], 2);
        denom += 2*p[i];
    }
    return num/denom;
}

SAXSDataset ScatteringHistogram::calc_guinier_approx() const {
    double Rg2 = calc_guinier_gyration_ratio_squared();

    const Axis& debye_axis = setting::axes::scattering_intensity_plot_axis;
    vector<double> Iq(debye_axis.bins, 0);
    double log_conv = log10(exp(1)); // we have to convert natural log to log10
    for (unsigned int i = 0; i < debye_axis.bins; i++) { // iterate through all q values
        Iq[i] = -pow(q[i], 2)*Rg2/3*log_conv;
    }

    return SAXSDataset(q, Iq, "q", "logI");
}

ScatteringHistogram& ScatteringHistogram::operator=(const ScatteringHistogram& h) {
    return operator=(ScatteringHistogram(h));
}

ScatteringHistogram& ScatteringHistogram::operator=(ScatteringHistogram&& h) {
    p = std::move(h.p);
    p_pp = std::move(h.p_pp);
    p_hh = std::move(h.p_hh);
    p_hp = std::move(h.p_hp);
    q = std::move(h.q);
    d = std::move(h.d);
    sinqd_table = std::move(h.sinqd_table);
    return *this;
}