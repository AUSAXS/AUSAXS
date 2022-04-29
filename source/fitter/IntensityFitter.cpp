#include <fitter/Fit.h>
#include <fitter/IntensityFitter.h>
#include <math/SimpleLeastSquares.h>
#include <math/CubicSpline.h>
#include <histogram/ScatteringHistogram.h>
#include <utility/Exceptions.h>

#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>

using std::string, std::vector, std::shared_ptr, std::unique_ptr;

shared_ptr<Fit> IntensityFitter::fit() {
    ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    auto f = std::bind(&IntensityFitter::chi2, this, std::placeholders::_1);
    ROOT::Math::Functor functor(f, 1); // declare the function to be minimized and its number of parameters
    minimizer->SetFunction(functor);
    minimizer->SetPrintLevel(0);
    minimizer->SetLimitedVariable(0, "c", 5, 1e-4, 0, 100); // scaling factor
    minimizer->Minimize();
    const double* res = minimizer->X();
    const double* err = minimizer->Errors();

    // apply c
    h.apply_water_scaling_factor(res[0]);
    vector<double> ym = h.calc_debye_scattering_intensity().get("I");
    vector<double> Im = splice(ym);

    // fit a, b
    SimpleLeastSquares fitter(Im, Io, sigma);
    std::shared_ptr<Fit> ab_fit = fitter.fit();

    // update fitter object
    bool converged = !minimizer->Status();
    std::map<string, double> pars = {{"c", res[0]}, {"a", ab_fit->params["a"]}, {"b", ab_fit->params["b"]}};
    std::map<string, double> errs = {{"c", err[0]}, {"a", ab_fit->errors["a"]}, {"b", ab_fit->errors["b"]}};
    double funcalls = minimizer->NCalls();
    fitted = std::make_shared<Fit>(pars, errs, chi2(res), qo.size()-2, funcalls, converged);

    minimizer->SetPrintLevel(3);
    minimizer->PrintResults();
    return fitted;
}

vector<shared_ptr<TGraph>> IntensityFitter::plot() {
    if (fitted == nullptr) {throw except::bad_order("Error in IntensityFitter::plot: Cannot plot before a fit has been made!");}

    double a = fitted->params["a"];
    double b = fitted->params["b"];
    double c = fitted->params["c"];

    h.apply_water_scaling_factor(c);
    vector<double> ym = h.calc_debye_scattering_intensity().get("I");
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

unique_ptr<TGraphErrors> IntensityFitter::plot_residuals() {
    if (fitted == nullptr) {throw except::bad_order("Error in IntensityFitter::plot_residuals: Cannot plot before a fit has been made!");}
 
    double a = fitted->params["a"];
    double b = fitted->params["b"];
    double c = fitted->params["c"];

    h.apply_water_scaling_factor(c);
    vector<double> ym = h.calc_debye_scattering_intensity().get("I");
    vector<double> Im = splice(ym);

    // calculate the residuals
    vector<double> residuals(qo.size());
    for (size_t i = 0; i < qo.size(); ++i) {
        residuals[i] = ((Io[i] - (a*Im[i]+b))/sigma[i]);
    }

    // prepare the TGraph
    vector<double> xerr(sigma.size(), 0);
    // unique_ptr<TGraphErrors> graph = std::make_unique<TGraphErrors>(qo.size(), &qo[0], &residuals[0], &xerr[0], &sigma[0]);
    unique_ptr<TGraphErrors> graph = std::make_unique<TGraphErrors>(qo.size(), &qo[0], &residuals[0], &xerr[0], &xerr[0]);
    return graph;
}

double IntensityFitter::chi2(const double* params) {
    double c = params[0];

    // apply c
    h.apply_water_scaling_factor(c);
    vector<double> ym = h.calc_debye_scattering_intensity().get("I");
    vector<double> Im = splice(ym);

    // fit a, b
    SimpleLeastSquares fitter(Im, Io, sigma);
    auto[a, b] = fitter.fit_params_only();

    // calculate chi2
    double chi = 0;
    for (size_t i = 0; i < qo.size(); i++) {
        double v = (Io[i] - (a*Im[i]+b))/sigma[i];
        chi += v*v;
    }
    return chi;
}

double IntensityFitter::get_intercept() {
    if (fitted == nullptr) {throw except::bad_order("Error in IntensityFitter::get_intercept: Cannot determine model intercept before a fit has been made!");}
 
    double a = fitted->params["a"];
    double b = fitted->params["b"];
    double c = fitted->params["c"];

    h.apply_water_scaling_factor(c);
    vector<double> ym = h.calc_debye_scattering_intensity().get("I");
    CubicSpline s(h.q, ym);
    return a*s.spline(0) + b;
}

SAXSDataset IntensityFitter::get_model_dataset() {
    if (fitted == nullptr) {throw except::bad_order("Error in IntensityFitter::get_model_dataset: Cannot determine model intercept before a fit has been made!");}
 
    double a = fitted->params["a"];
    double b = fitted->params["b"];
    double c = fitted->params["c"];

    h.apply_water_scaling_factor(c);
    vector<double> ym = h.calc_debye_scattering_intensity().get("I");
    vector<double> Im = splice(ym);
    std::transform(Im.begin(), Im.end(), Im.begin(), [&a, &b] (double I) {return I*a+b;});

    return SAXSDataset(qo, Im, "q", "I"); 
}

SAXSDataset IntensityFitter::get_model_dataset(const vector<double>& q) {
    if (fitted == nullptr) {throw except::bad_order("Error in IntensityFitter::get_model_dataset: Cannot determine model intercept before a fit has been made!");}
 
    double a = fitted->params["a"];
    double b = fitted->params["b"];
    double c = fitted->params["c"];

    h.apply_water_scaling_factor(c);
    SAXSDataset data = h.calc_debye_scattering_intensity(q);
    std::transform(data.y.begin(), data.y.end(), data.y.begin(), [&a, &b] (double I) {return I*a+b;});
    return data;
}

SAXSDataset IntensityFitter::get_dataset() const {
    SAXSDataset data(qo, Io, sigma);
    data.xlabel = "q";
    data.ylabel = "I";
    return data;
}