#include <iostream>
#include <fstream>
#include <stdexcept>
#include <tuple>
#include <map>

#include "fitter/SimpleIntensityFitter.h"
#include "math/CubicSpline.h"
#include "settings.h"
#include "ScatteringHistogram.h"
#include "Exceptions.h"

#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <TGraph.h>
#include <TGraphErrors.h>

using std::string, std::vector, std::shared_ptr, std::unique_ptr;

SimpleIntensityFitter::SimpleIntensityFitter(const ScatteringHistogram& data, const ScatteringHistogram& model) : h(data) {
    model_setup(model);
}

SimpleIntensityFitter::SimpleIntensityFitter(const ScatteringHistogram& model) {
    model_setup(model);
}

void SimpleIntensityFitter::model_setup(const ScatteringHistogram& model) {
    double qmin = setting::fit::q_low;
    double qmax = setting::fit::q_high;

    // calculate q & sigma
    qo.resize(100);
    sigma.resize(100); 
    double width = (qmax-qmin)/100;
    for (unsigned int i = 0; i < 100; i++) {
        qo[i] = qmin + i*width;
        sigma[i] = 1;
    }

    // calculate the intensity
    Io = model.calc_debye_scattering_intensity(qo);
}

shared_ptr<Fitter::Fit> SimpleIntensityFitter::fit() {
    vector<double> ym = h.calc_debye_scattering_intensity();
    vector<double> Im = splice(ym);

    SimpleLeastSquares fitter(Im, Io, sigma);
    fitted = fitter.fit();
    return fitted;
}

vector<shared_ptr<TGraph>> SimpleIntensityFitter::plot() {
    if (fitted == nullptr) {throw except::bad_order("Error in IntensityFitter::plot: Cannot plot before a fit has been made!");}

    double a = fitted->params["a"];
    double b = fitted->params["b"];

    vector<double> ym = h.calc_debye_scattering_intensity();
    vector<double> Im = splice(ym);

    // calculate the scaled I model values
    vector<double> I_scaled(qo.size()); // spliced data
    vector<double> ym_scaled(ym.size()); // original scaled data
    std::transform(Im.begin(), Im.end(), I_scaled.begin(), [&a, &b] (double I) {return I*a+b;});
    std::transform(ym.begin(), ym.end(), ym_scaled.begin(), [&a, &b] (double I) {return I*a+b;});

    // prepare the TGraphs
    vector<double> xerr(sigma.size(), 0);
    vector<shared_ptr<TGraph>> graphs(3);
    graphs[0] = std::make_shared<TGraph>(qo.size(), &qo[0], &I_scaled[0]);
    graphs[1] = std::make_shared<TGraph>(h.q.size(), &h.q[0], &ym_scaled[0]);
    graphs[2] = std::make_shared<TGraphErrors>(qo.size(), &qo[0], &Io[0], &xerr[0], &sigma[0]);
    return graphs;
}

unique_ptr<TGraphErrors> SimpleIntensityFitter::plot_residuals() {
    if (fitted == nullptr) {throw except::bad_order("Error in IntensityFitter::plot_residuals: Cannot plot before a fit has been made!");}
 
    double a = fitted->params["a"];
    double b = fitted->params["b"];

    vector<double> ym = h.calc_debye_scattering_intensity();
    vector<double> Im = splice(ym);

    // calculate the residuals
    vector<double> residuals(qo.size());
    for (size_t i = 0; i < qo.size(); ++i) {
        residuals[i] = ((Io[i] - a*Im[i]-b)/sigma[i]);
    }

    // prepare the TGraph
    vector<double> xerr(sigma.size(), 0);
    unique_ptr<TGraphErrors> graph = std::make_unique<TGraphErrors>(qo.size(), &qo[0], &residuals[0], &xerr[0], &sigma[0]);
    return graph;
}

void SimpleIntensityFitter::set_scattering_hist(ScatteringHistogram&& h) {
    this->h = std::move(h);
}

void SimpleIntensityFitter::set_scattering_hist(const ScatteringHistogram& h) {
    this->h = h;
}

double SimpleIntensityFitter::chi2(const double* params) {
    vector<double> ym = h.calc_debye_scattering_intensity();
    vector<double> Im = splice(ym);

    // fit a, b
    SimpleLeastSquares fitter(Im, Io, sigma);
    auto[a, b] = fitter.fit_params_only();

    // calculate chi2
    double chi = 0;
    for (size_t i = 0; i < qo.size(); i++) {
        chi += pow((Io[i] - a*Im[i]-b)/sigma[i], 2);
    }
    return chi;
}

void SimpleIntensityFitter::setup(string file) {
    std::tie(qo, Io, sigma) = read(file); // read observed values from input file
}

vector<double> SimpleIntensityFitter::splice(const vector<double>& ym) const {
    vector<double> Im = vector<double>(qo.size()); // spliced model values
    CubicSpline s(h.q, ym);
    for (size_t i = 0; i < qo.size(); ++i) {
        Im[i] = s.spline(qo[i]);
    }
    return Im;
}

std::tuple<vector<double>, vector<double>, vector<double>> SimpleIntensityFitter::read(string file) const {
    // check if file was succesfully opened
    std::ifstream input(file);
    if (!input.is_open()) {throw std::ios_base::failure("Error in IntensityFitter::read: Could not open file \"" + file + "\"");}

    vector<double> q;
    vector<double> I;
    vector<double> sigma;
    string line; // placeholder for the current line
    while(getline(input, line)) {
        if (line[0] == ' ') {line = line.substr(1);} // fix leading space
        vector<string> tokens;
        boost::split(tokens, line, boost::is_any_of(" ,\t")); // spaces, commas, and tabs can all be used as separators (but not a mix of them)

        // determine if we are in some sort of header
        if (tokens.size() < 3 || tokens.size() > 4) {continue;} // too many separators
        bool skip = false;
        for (int i = 0; i < 3; i++) { // check if they are numbers
            if (!tokens[i].empty() && tokens[i].find_first_not_of("0123456789-.Ee") != string::npos) {skip = true;}
        }
        if (skip) {continue;}

        // now we are most likely beyond any headers
        double _q, _I, _sigma;
        _q = std::stod(tokens[0]); // we know for sure that the strings are convertible to numbers (boost check)
        _I = std::stod(tokens[1]);
        _sigma = std::stod(tokens[2]);

        if (_q > 10) {continue;} // probably not a q-value if it's larger than 10

        // check user-defined limits
        if (_q < setting::fit::q_low) {continue;}
        if (_q > setting::fit::q_high) {continue;}

        // add the values to our vectors
        q.push_back(_q);
        I.push_back(_I);
        sigma.push_back(_sigma); 
    }
    return std::make_tuple(q, I, sigma);
}