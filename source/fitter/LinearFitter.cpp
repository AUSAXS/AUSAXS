#include <iostream>
#include <fstream>

#include <fitter/LinearFitter.h>
#include <math/CubicSpline.h>
#include <math/SimpleLeastSquares.h>
#include <hist/ScatteringHistogram.h>
#include <utility/Exceptions.h>
#include <settings/HistogramSettings.h>

using namespace fitter;

LinearFitter::LinearFitter(std::string input) {setup(input);}
LinearFitter::LinearFitter(std::string input, const hist::ScatteringHistogram& h) : h(h) {setup(input);}
LinearFitter::LinearFitter(std::string input, hist::ScatteringHistogram&& h) : h(std::move(h)) {setup(input);}
LinearFitter::LinearFitter(const SimpleDataset& data) : data(data) {}
LinearFitter::LinearFitter(const SimpleDataset& data, const hist::ScatteringHistogram& hist) : data(data), h(hist) {}
LinearFitter::LinearFitter(const hist::ScatteringHistogram& data, const hist::ScatteringHistogram& model) : LinearFitter(data, model, Limit(settings::axes::qmin, settings::axes::qmax)) {}
LinearFitter::LinearFitter(const hist::ScatteringHistogram& model) : LinearFitter(model, Limit(settings::axes::qmin, settings::axes::qmax)) {}
LinearFitter::LinearFitter(const hist::ScatteringHistogram& data, const hist::ScatteringHistogram& model, const Limit& limits) : h(data) {
    model_setup(model, limits);
}
LinearFitter::LinearFitter(const hist::ScatteringHistogram& model, const Limit& limits) {
    model_setup(model, limits);
}

void LinearFitter::model_setup(const hist::ScatteringHistogram& model, const Limit& limits) {
    throw except::not_implemented("LinearFitter::model_setup: Not implemented yet!");
    // data = model.calc_debye_scattering_intensity();
    // data.reduce(settings::fit::N, true);
    // data.limit_x(limits);
    // data.simulate_errors();
    // if (I0 > 0) {data.normalize(I0);}
    // if (settings::em::simulation::noise) {data.simulate_noise();}
}

double LinearFitter::fit_only() {
    return chi2({});
}

std::shared_ptr<Fit> LinearFitter::fit() {
    std::vector<double> ym = h.calc_debye_scattering_intensity().col("I");
    std::vector<double> Im = splice(ym);

    // we want to fit a*Im + b to Io
    SimpleDataset fit_data(Im, data.y(), data.yerr());
    if (I0 > 0) {fit_data.normalize(I0);}

    SimpleLeastSquares fitter(fit_data);
    fitted = fitter.fit();
    fitted->add_plots(*this);

    return fitted;
}

void LinearFitter::normalize_intensity(double new_I0) {
    if (I0 < 0) {data.normalize(new_I0);} // if y0 has not been set yet, we must rescale the data
    I0 = new_I0;
}

Fit::Plots LinearFitter::plot() {
    if (fitted == nullptr) {throw except::bad_order("IntensityFitter::plot: Cannot plot before a fit has been made!");}

    double a = fitted->get_parameter("a").value;
    double b = fitted->get_parameter("b").value;

    std::vector<double> ym = h.calc_debye_scattering_intensity().col("I");
    std::vector<double> Im = splice(ym);

    // if we have a I0, we need to rescale the data
    // double factor = I0/ym[0];
    // std::transform(Im.begin(), Im.end(), Im.begin(), [&factor] (double y) {return factor*y;});

    // calculate the scaled I model values
    std::vector<double> I_scaled(data.size()); // spliced data
    std::vector<double> ym_scaled(ym.size()); // original scaled data
    std::transform(Im.begin(), Im.end(), I_scaled.begin(), [&a, &b] (double I) {return I*a+b;});
    std::transform(ym.begin(), ym.end(), ym_scaled.begin(), [&a, &b] (double I) {return I*a+b;});

    // prepare the TGraphs
    Fit::Plots graphs;
    graphs.intensity_interpolated = SimpleDataset(data.x(), I_scaled);
    graphs.intensity = SimpleDataset(h.q, ym_scaled);
    graphs.data = SimpleDataset(data.x(), data.y(), data.yerr());

    auto lim = graphs.data.get_xlimits();
    lim.expand(0.05);
    graphs.intensity.limit_x(lim);
    return graphs;
}

SimpleDataset LinearFitter::plot_residuals() {
    if (fitted == nullptr) {throw except::bad_order("IntensityFitter::plot_residuals: Cannot plot before a fit has been made!");}
 
    double a = fitted->get_parameter("a").value;
    double b = fitted->get_parameter("b").value;

    std::vector<double> ym = h.calc_debye_scattering_intensity().col("I");
    std::vector<double> Im = splice(ym);

    // calculate the residuals
    std::vector<double> residuals(data.size());
    for (size_t i = 0; i < data.size(); ++i) {
        residuals[i] = ((data.y(i) - a*Im[i]-b)/data.yerr(i));
    }

    // prepare the dataset
    std::vector<double> xerr(data.size(), 0);
    return Dataset2D(data.x(), residuals, xerr, data.yerr());
}

void LinearFitter::set_scattering_hist(hist::ScatteringHistogram&& h) {
    this->h = std::move(h);
}

void LinearFitter::set_scattering_hist(const hist::ScatteringHistogram& h) {
    this->h = h;
}

double LinearFitter::chi2(const std::vector<double>&) {
    std::vector<double> ym = h.calc_debye_scattering_intensity().col("I");
    std::vector<double> Im = splice(ym);

    // we want to fit a*Im + b to Io
    SimpleDataset fit_data(Im, data.y(), data.yerr());
    if (I0 > 0) {fit_data.normalize(I0);}

    SimpleLeastSquares fitter(fit_data);
    return fitter.fit_only();
}

void LinearFitter::setup(std::string file) {
    data = SimpleDataset(file); // read observed values from input file
}

std::vector<double> LinearFitter::splice(const std::vector<double>& ym) const {
    std::vector<double> Im(data.size()); // spliced model values
    CubicSpline s(h.q, ym);
    for (size_t i = 0; i < data.size(); ++i) {
        Im[i] = s.spline(data.x(i));
    }
    return Im;
}

unsigned int LinearFitter::degrees_of_freedom() const {
    return data.size() - 2;
}

unsigned int LinearFitter::dof() const {
    return degrees_of_freedom();
}

std::shared_ptr<Fit> LinearFitter::get_fit() const {
    if (fitted == nullptr) {throw except::bad_order("LinearFitter::get_fit: Cannot get the fit results before a fit has been made!");}
    return fitted;
}