#include <fitter/Fit.h>
#include <fitter/IntensityFitter.h>
#include <math/SimpleLeastSquares.h>
#include <math/CubicSpline.h>
#include <hist/ScatteringHistogram.h>
#include <utility/Exceptions.h>
#include <mini/all.h>
#include <plots/all.h>

IntensityFitter::IntensityFitter(const hist::ScatteringHistogram& model, const Limit& limits) : SimpleIntensityFitter(model, limits) {}

std::shared_ptr<Fit> IntensityFitter::fit() {
    std::function<double(std::vector<double>)> f = std::bind(&IntensityFitter::chi2, this, std::placeholders::_1);
    auto mini = mini::create_minimizer(fit_type, f, guess, setting::em::evals);
    auto res = mini->minimize();

    // mini::Scan mini(f, guess, 1000);
    // auto res = mini.minimize();
    // auto d = mini.get_evaluated_points();
    // d.add_plot_options("points", {{"xlabel", "c"}, {"ylabel", "chi2"}});
    // plots::PlotDataset::quick_plot(d, "scan.pdf");

    // apply c
    h.apply_water_scaling_factor(res.get_parameter("c").value);
    std::vector<double> ym = h.calc_debye_scattering_intensity().col("I");
    std::vector<double> Im = splice(ym);

    // we want to fit a*Im + b to Io
    SimpleDataset fit_data(Im, data.y(), data.yerr());
    if (I0 > 0) {fit_data.normalize(I0);}

    SimpleLeastSquares fitter(fit_data);
    std::shared_ptr<Fit> ab_fit = fitter.fit();

    // update fitter object
    fitted = std::make_shared<Fit>(res, res.fval, data.size()-1); // start with the fit performed here
    fitted->add_fit(ab_fit);                                      // add the a,b inner fit
    fitted->add_plots(*this);                                     // make the result plottable
    fitted->evaluated_points = mini->get_evaluated_points();      // add the evaluated points

    return fitted;
}

Fit::Plots IntensityFitter::plot() {
    if (fitted == nullptr) {throw except::bad_order("IntensityFitter::plot: Cannot plot before a fit has been made!");}

    double a = fitted->get_parameter("a").value;
    double b = fitted->get_parameter("b").value;
    double c = fitted->get_parameter("c").value;

    h.apply_water_scaling_factor(c);
    std::vector<double> ym = h.calc_debye_scattering_intensity().col("I");
    std::vector<double> Im = splice(ym);

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

SimpleDataset IntensityFitter::plot_residuals() {
    if (fitted == nullptr) {throw except::bad_order("IntensityFitter::plot_residuals: Cannot plot before a fit has been made!");}
 
    double a = fitted->get_parameter("a").value;
    double b = fitted->get_parameter("b").value;
    double c = fitted->get_parameter("c").value;

    h.apply_water_scaling_factor(c);
    std::vector<double> ym = h.calc_debye_scattering_intensity().col("I");
    std::vector<double> Im = splice(ym);

    // calculate the residuals
    std::vector<double> residuals(data.size());
    for (size_t i = 0; i < data.size(); ++i) {
        residuals[i] = ((data.y(i) - (a*Im[i]+b))/data.yerr(i));
    }

    // prepare the TGraph
    std::vector<double> xerr(data.size(), 0);
    return Dataset2D(data.x(), residuals, xerr, data.yerr());
}

double IntensityFitter::chi2(std::vector<double> params) {
    double c = params[0];

    // apply c
    h.apply_water_scaling_factor(c);
    std::vector<double> ym = h.calc_debye_scattering_intensity().col("I");
    std::vector<double> Im = splice(ym);

    // we want to fit a*Im + b to Io
    SimpleDataset fit_data(Im, data.y(), data.yerr());
    if (I0 > 0) {fit_data.normalize(I0);}

    SimpleLeastSquares fitter(fit_data);
    auto[a, b] = fitter.fit_params_only();

    // calculate chi2
    double chi = 0;
    for (size_t i = 0; i < data.size(); i++) {
        double v = (data.y(i) - (a*Im[i]+b))/data.yerr(i);
        chi += v*v;
    }

    return chi;
}

double IntensityFitter::get_intercept() {
    if (fitted == nullptr) {throw except::bad_order("IntensityFitter::get_intercept: Cannot determine model intercept before a fit has been made!");}
 
    double a = fitted->get_parameter("a").value;
    double b = fitted->get_parameter("b").value;
    double c = fitted->get_parameter("c").value;

    h.apply_water_scaling_factor(c);
    std::vector<double> ym = h.calc_debye_scattering_intensity().col("I");
    CubicSpline s(h.q, ym);
    return a*s.spline(0) + b;
}

SimpleDataset IntensityFitter::get_model_dataset() {
    if (fitted == nullptr) {throw except::bad_order("IntensityFitter::get_model_dataset: Cannot determine model intercept before a fit has been made!");}
 
    double a = fitted->get_parameter("a").value;
    double b = fitted->get_parameter("b").value;
    double c = fitted->get_parameter("c").value;

    h.apply_water_scaling_factor(c);
    std::vector<double> ym = h.calc_debye_scattering_intensity().col("I");
    std::vector<double> Im = splice(ym);
    std::transform(Im.begin(), Im.end(), Im.begin(), [&a, &b] (double I) {return I*a+b;});

    return SimpleDataset(data.x(), Im, "q", "I"); 
}

SimpleDataset IntensityFitter::get_model_dataset(const std::vector<double>& q) {
    if (fitted == nullptr) {throw except::bad_order("IntensityFitter::get_model_dataset: Cannot determine model intercept before a fit has been made!");}
 
    double a = fitted->get_parameter("a").value;
    double b = fitted->get_parameter("b").value;
    double c = fitted->get_parameter("c").value;

    h.apply_water_scaling_factor(c);
    SimpleDataset model = h.calc_debye_scattering_intensity(q);
    auto y = model.y();
    std::transform(y.begin(), y.end(), y.begin(), [&a, &b] (double I) {return I*a+b;});
    return model;
}

SimpleDataset IntensityFitter::get_dataset() const {
    return data;
}

void IntensityFitter::set_guess(mini::Parameter guess) {
    this->guess = guess;
}