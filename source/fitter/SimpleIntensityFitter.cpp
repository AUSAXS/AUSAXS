#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>

#include <fitter/SimpleIntensityFitter.h>
#include <math/CubicSpline.h>
#include <math/SimpleLeastSquares.h>
#include <utility/Settings.h>
#include <histogram/ScatteringHistogram.h>
#include <utility/Exceptions.h>

#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>

using std::string, std::vector, std::shared_ptr, std::unique_ptr;

SimpleIntensityFitter::SimpleIntensityFitter(const hist::ScatteringHistogram& data, const hist::ScatteringHistogram& model, const Limit& limits) : h(data) {
    model_setup(model, limits);
}

SimpleIntensityFitter::SimpleIntensityFitter(const hist::ScatteringHistogram& model, const Limit& limits) {
    model_setup(model, limits);
}

void SimpleIntensityFitter::model_setup(const hist::ScatteringHistogram& model, const Limit& limits) {
    data = model.calc_debye_scattering_intensity();
    data.reduce(setting::fit::N, true);
    data.limit(limits);
    data.simulate_errors();
    if (I0 > 0) {data.normalize(I0);}
    if (setting::em::simulation::noise) {data.simulate_noise();}
}

shared_ptr<Fit> SimpleIntensityFitter::fit() {
    vector<double> ym = h.calc_debye_scattering_intensity().get("I");
    vector<double> Im = splice(ym);

    // we want to fit a*Im + b to Io
    Dataset fit_data(Im, data.y, data.yerr);
    if (I0 > 0) {fit_data.normalize(I0);}

    SimpleLeastSquares fitter(fit_data);
    fitted = fitter.fit();

    return fitted;
}

void SimpleIntensityFitter::normalize_intensity(double new_I0) {
    if (I0 < 0) {data.normalize(new_I0);} // if y0 has not been set yet, we must rescale the data
    I0 = new_I0;
}

Fit::Plots SimpleIntensityFitter::plot() {
    if (fitted == nullptr) {throw except::bad_order("Error in IntensityFitter::plot: Cannot plot before a fit has been made!");}

    double a = fitted->get_parameter("a").value;
    double b = fitted->get_parameter("b").value;

    vector<double> ym = h.calc_debye_scattering_intensity().get("I");
    vector<double> Im = splice(ym);

    // if we have a I0, we need to rescale the data
    // double factor = I0/ym[0];
    // std::transform(Im.begin(), Im.end(), Im.begin(), [&factor] (double y) {return factor*y;});

    // calculate the scaled I model values
    vector<double> I_scaled(data.size()); // spliced data
    vector<double> ym_scaled(ym.size()); // original scaled data
    std::transform(Im.begin(), Im.end(), I_scaled.begin(), [&a, &b] (double I) {return I*a+b;});
    std::transform(ym.begin(), ym.end(), ym_scaled.begin(), [&a, &b] (double I) {return I*a+b;});

    // prepare the TGraphs
    vector<double> xerr(data.size(), 0);
    Fit::Plots graphs;
    graphs.intensity_interpolated = SAXSDataset(data.x, I_scaled);
    graphs.intensity = SAXSDataset(h.q, ym_scaled);
    graphs.data = SAXSDataset(data.x, data.y, xerr, data.yerr);
    return graphs;
}

Dataset SimpleIntensityFitter::plot_residuals() {
    if (fitted == nullptr) {throw except::bad_order("Error in IntensityFitter::plot_residuals: Cannot plot before a fit has been made!");}
 
    double a = fitted->get_parameter("a").value;
    double b = fitted->get_parameter("b").value;

    vector<double> ym = h.calc_debye_scattering_intensity().get("I");
    vector<double> Im = splice(ym);

    // calculate the residuals
    vector<double> residuals(data.size());
    for (size_t i = 0; i < data.size(); ++i) {
        residuals[i] = ((data.y[i] - a*Im[i]-b)/data.yerr[i]);
    }

    // prepare the dataset
    vector<double> xerr(data.size(), 0);
    return Dataset(data.x, residuals, xerr, data.yerr);
}

void SimpleIntensityFitter::set_scattering_hist(hist::ScatteringHistogram&& h) {
    this->h = std::move(h);
}

void SimpleIntensityFitter::set_scattering_hist(const hist::ScatteringHistogram& h) {
    this->h = h;
}

double SimpleIntensityFitter::chi2(const double*) {
    throw except::invalid_operation("Error in SimpleIntensityFitter::chi2: Not implemented.");
    // vector<double> ym = h.calc_debye_scattering_intensity().get("I");
    // vector<double> Im = splice(ym);

    // // fit a, b
    // SimpleLeastSquares fitter(Im, Io, sigma);
    // auto[a, b] = fitter.fit_params_only();

    // // calculate chi2
    // double chi = 0;
    // for (size_t i = 0; i < qo.size(); i++) {
    //     chi += pow((Io[i] - a*Im[i]-b)/sigma[i], 2);
    // }
    // return chi;
}

void SimpleIntensityFitter::setup(string file) {
    data = read(file); // read observed values from input file
}

vector<double> SimpleIntensityFitter::splice(const vector<double>& ym) const {
    vector<double> Im = vector<double>(data.size()); // spliced model values
    CubicSpline s(h.q, ym);
    for (size_t i = 0; i < data.size(); ++i) {
        Im[i] = s.spline(data.x[i]);
    }
    return Im;
}

SAXSDataset SimpleIntensityFitter::read(string file) const {
    // check if file was succesfully opened
    std::ifstream input(file);
    if (!input.is_open()) {throw std::ios_base::failure("Error in IntensityFitter::read: Could not open file \"" + file + "\"");}

    SAXSDataset temp;
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
        temp.x.push_back(_q);
        temp.y.push_back(_I);
        temp.yerr.push_back(_sigma); 
    }
    return temp;
}

unsigned int SimpleIntensityFitter::degrees_of_freedom() const {
    return data.size() - 2;
}

unsigned int SimpleIntensityFitter::dof() const {
    return degrees_of_freedom();
}

std::shared_ptr<Fit> SimpleIntensityFitter::get_fit() const {
    if (fitted == nullptr) {throw except::bad_order("Error in SimpleIntensityFitter::get_fit: Cannot get the fit results before a fit has been made!");}
    return fitted;
}