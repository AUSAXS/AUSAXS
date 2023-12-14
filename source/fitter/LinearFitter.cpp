#include <fitter/LinearFitter.h>
#include <fitter/FitPlots.h>
#include <math/CubicSpline.h>
#include <math/SimpleLeastSquares.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <hist/Histogram.h>
#include <utility/Exceptions.h>
#include <settings/HistogramSettings.h>
#include <mini/detail/FittedParameter.h>
#include <mini/detail/Evaluation.h>
#include <dataset/Dataset2D.h>
#include <settings/EMSettings.h>

#include <iostream>
#include <fstream>

using namespace fitter;

LinearFitter::LinearFitter(LinearFitter&& other) : fitted(std::move(other.fitted)), data(std::move(other.data)), I0(other.I0), h(std::move(other.h)) {}
LinearFitter::LinearFitter(const io::ExistingFile& input) {setup(input);}
LinearFitter::LinearFitter(const io::ExistingFile& input, std::unique_ptr<hist::DistanceHistogram> h) : h(std::move(h)) {setup(input);}
LinearFitter::LinearFitter(const SimpleDataset& data) : data(data) {}
LinearFitter::LinearFitter(const SimpleDataset& data, std::unique_ptr<hist::DistanceHistogram> h) : data(data), h(std::move(h)) {}
LinearFitter::LinearFitter(std::unique_ptr<hist::DistanceHistogram> data, std::unique_ptr<hist::DistanceHistogram> model) : LinearFitter(std::move(data), std::move(model), Limit(settings::axes::qmin, settings::axes::qmax)) {}
LinearFitter::LinearFitter(std::unique_ptr<hist::DistanceHistogram> data, std::unique_ptr<hist::DistanceHistogram> model, const Limit& limits) : h(std::move(data)) {
    model_setup(std::move(model), limits);
}
LinearFitter::LinearFitter(std::unique_ptr<hist::DistanceHistogram> model) : LinearFitter(std::move(model), Limit(settings::axes::qmin, settings::axes::qmax)) {}
LinearFitter::LinearFitter(std::unique_ptr<hist::DistanceHistogram> model, const Limit& limits) {
    model_setup(std::move(model), limits);
}

void LinearFitter::model_setup(std::unique_ptr<hist::DistanceHistogram> model, const Limit& limits) {
    data = model->debye_transform();
    data.limit_x(limits);
    data.simulate_errors();
    if (I0 > 0) {data.normalize(I0);}
    if (settings::em::simulation::noise) {data.simulate_noise();}
}

double LinearFitter::fit_chi2_only() {
    return chi2({});
}

std::shared_ptr<Fit> LinearFitter::fit() {
    std::vector<double> ym = h->debye_transform().get_counts();
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

FitPlots LinearFitter::plot() {
    if (fitted == nullptr) {throw except::bad_order("IntensityFitter::plot: Cannot plot before a fit has been made!");}

    double a = fitted->get_parameter("a").value;
    double b = fitted->get_parameter("b").value;

    std::vector<double> ym = h->debye_transform().get_counts();
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
    FitPlots graphs;
    graphs.intensity_interpolated = SimpleDataset(data.x(), I_scaled);
    graphs.intensity = SimpleDataset(h->get_q_axis(), ym_scaled);
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

    std::vector<double> ym = h->debye_transform().get_counts();
    std::vector<double> Im = splice(ym);

    // calculate the residuals
    std::vector<double> residuals(data.size());
    for (int i = 0; i < static_cast<int>(data.size()); ++i) {
        residuals[i] = ((data.y(i) - a*Im[i]-b)/data.yerr(i));
    }

    // prepare the dataset
    std::vector<double> xerr(data.size(), 0);
    return Dataset2D(data.x(), residuals, xerr, data.yerr());
}

void LinearFitter::set_scattering_hist(std::unique_ptr<hist::DistanceHistogram> h) {
    this->h = std::move(h);
}

observer_ptr<hist::DistanceHistogram> LinearFitter::get_scattering_hist() {
    return h.get();
}

double LinearFitter::chi2(const std::vector<double>&) {
    std::vector<double> ym = h->debye_transform().get_counts();
    std::vector<double> Im = splice(ym);

    // we want to fit a*Im + b to Io
    SimpleDataset fit_data(Im, data.y(), data.yerr());
    if (I0 > 0) {fit_data.normalize(I0);}

    SimpleLeastSquares fitter(fit_data);
    return fitter.fit_chi2_only();
}

void LinearFitter::setup(const io::ExistingFile& file) {
    data = SimpleDataset(file); // read observed values from input file
}

std::vector<double> LinearFitter::splice(const std::vector<double>& ym) const {
    std::vector<double> Im(data.size()); // spliced model values
    math::CubicSpline s(h->get_q_axis(), ym);
    for (unsigned int i = 0; i < data.size(); ++i) {
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

void LinearFitter::operator=(LinearFitter&& other) {
    fitted = std::move(other.fitted);
    data = std::move(other.data);
    I0 = other.I0;
    h = std::move(other.h);
}