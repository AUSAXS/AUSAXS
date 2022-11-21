#include <iostream>
#include <fstream>

#include <fitter/SimpleIntensityFitter.h>
#include <math/CubicSpline.h>
#include <math/SimpleLeastSquares.h>
#include <utility/Settings.h>
#include <hist/ScatteringHistogram.h>
#include <utility/Exceptions.h>

SimpleIntensityFitter::SimpleIntensityFitter(const hist::ScatteringHistogram& data, const hist::ScatteringHistogram& model, const Limit& limits) : h(data) {
    model_setup(model, limits);
}

SimpleIntensityFitter::SimpleIntensityFitter(const hist::ScatteringHistogram& model, const Limit& limits) {
    model_setup(model, limits);
}

void SimpleIntensityFitter::model_setup(const hist::ScatteringHistogram& model, const Limit& limits) {
    data = model.calc_debye_scattering_intensity();
    data.reduce(setting::fit::N, true);
    data.limit_x(limits);
    data.simulate_errors();
    if (I0 > 0) {data.normalize(I0);}
    if (setting::em::simulation::noise) {data.simulate_noise();}
}

std::shared_ptr<Fit> SimpleIntensityFitter::fit() {
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

void SimpleIntensityFitter::normalize_intensity(double new_I0) {
    if (I0 < 0) {data.normalize(new_I0);} // if y0 has not been set yet, we must rescale the data
    I0 = new_I0;
}

Fit::Plots SimpleIntensityFitter::plot() {
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

SimpleDataset SimpleIntensityFitter::plot_residuals() {
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

void SimpleIntensityFitter::set_scattering_hist(hist::ScatteringHistogram&& h) {
    this->h = std::move(h);
}

void SimpleIntensityFitter::set_scattering_hist(const hist::ScatteringHistogram& h) {
    this->h = h;
}

double SimpleIntensityFitter::chi2(std::vector<double>) {
    throw except::invalid_operation("SimpleIntensityFitter::chi2: Not implemented.");
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

void SimpleIntensityFitter::setup(std::string file) {
    data = SimpleDataset(file); // read observed values from input file
}

std::vector<double> SimpleIntensityFitter::splice(const std::vector<double>& ym) const {
    std::vector<double> Im(data.size()); // spliced model values
    CubicSpline s(h.q, ym);
    for (size_t i = 0; i < data.size(); ++i) {
        Im[i] = s.spline(data.x(i));
    }
    return Im;
}

unsigned int SimpleIntensityFitter::degrees_of_freedom() const {
    return data.size() - 2;
}

unsigned int SimpleIntensityFitter::dof() const {
    return degrees_of_freedom();
}

std::shared_ptr<Fit> SimpleIntensityFitter::get_fit() const {
    if (fitted == nullptr) {throw except::bad_order("SimpleIntensityFitter::get_fit: Cannot get the fit results before a fit has been made!");}
    return fitted;
}