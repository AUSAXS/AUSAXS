#include "Fitter.h"
#include "math/CubicSpline.h"
#include "settings.h"
#include "SimpleLeastSquares.h"

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <tuple>
#include <map>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>

#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <TGraph.h>
#include <TGraphErrors.h>

using std::string, std::vector, std::shared_ptr, std::unique_ptr;

struct bad_order_except : public std::exception {
    bad_order_except(const char* msg) : msg(msg) {}
    const char* what() const throw() {return msg;}
    const char* msg;
};

class IntensityFitter : public Fitter {
public: 
    /**
     * @brief Prepare a fit of the measured values in @a input to the model described by @a q and @a I.
     * @param input the path to the file containing the measured values. 
     * @param q the model q values.
     * @param I the model I values. 
     */
    IntensityFitter(string input, vector<double>& q, vector<double>& I) : xm(q), ym(I) {setup(input, q, I);}
    ~IntensityFitter() override {}

    /**
     * @brief Perform the fit.
     * @return A Fit object containing various information about the fit. Note that the fitted scaling parameter is a = c/M*r_e^2 and b = background
     */
    // shared_ptr<Fit> fit() override {
    //     ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2");
    //     auto f = std::bind(&IntensityFitter::chi2, this, std::placeholders::_1);
    //     ROOT::Math::Functor functor(f, 2); // declare the function to be minimized and its number of parameters
    //     minimizer->SetFunction(functor);
    //     minimizer->SetVariable(0, "a", 0, 1e-4); // scaling factor
    //     minimizer->SetVariable(1, "b", 0, 1e-4); // background
    //     minimizer->Minimize();
    //     const double* res = minimizer->X();
    //     const double* err = minimizer->Errors();

    //     bool converged = !minimizer->Status();
    //     std::map<string, double> pars = {{"a", res[0]}, {"b", res[1]}};
    //     std::map<string, double> errs = {{"a", err[0]}, {"b", err[1]}};
    //     double funcalls = minimizer->NCalls();
    //     fitted = std::make_shared<Fit>(pars, errs, chi2(res), qo.size()-2, funcalls, converged);
    //     minimizer->PrintResults();
    //     return fitted;
    // }

    shared_ptr<Fit> fit() override {
        SimpleLeastSquares fitter(Im, Io, sigma);
        fitted = fitter.fit();
        return fitted;
    }

    vector<shared_ptr<TGraph>> plot() const {
        if (fitted == nullptr) {throw bad_order_except("Error in IntensityFitter::plot: Cannot plot before a fit has been made!");}

        // calculate the scaled I model values
        double a = fitted->params["a"];
        double b = fitted->params["b"];
        vector<double> I_scaled(qo.size()); // spliced scaled data
        vector<double> ym_scaled(ym.size()); // original scaled data
        std::transform(Im.begin(), Im.end(), I_scaled.begin(), [&a, &b] (double I) {return I*a+b;});
        std::transform(ym.begin(), ym.end(), ym_scaled.begin(), [&a, &b] (double I) {return I*a+b;});

        // prepare the TGraphs
        vector<double> xerr(sigma.size(), 0);
        vector<shared_ptr<TGraph>> graphs(3);
        graphs[0] = std::make_shared<TGraph>(qo.size(), &qo[0], &I_scaled[0]);
        graphs[1] = std::make_shared<TGraph>(xm.size(), &xm[0], &ym_scaled[0]);
        graphs[2] = std::make_shared<TGraphErrors>(qo.size(), &qo[0], &Io[0], &xerr[0], &sigma[0]);
        return graphs;
    }

    unique_ptr<TGraphErrors> plot_residuals() {
        if (fitted == nullptr) {throw bad_order_except("Error in IntensityFitter::plot_residuals: Cannot plot before a fit has been made!");}

        // calculate the residuals
        double a = fitted->params["a"];
        double b = fitted->params["b"];
        vector<double> residuals(qo.size());
        for (size_t i = 0; i < qo.size(); ++i) {
            residuals[i] = ((Io[i] - a*Im[i]-b)/sigma[i]);
        }

        // prepare the TGraph
        vector<double> xerr(sigma.size(), 0);
        unique_ptr<TGraphErrors> graph = std::make_unique<TGraphErrors>(qo.size(), &qo[0], &residuals[0], &xerr[0], &sigma[0]);
        return graph;
    }

private: 
    shared_ptr<Fit> fitted;
    vector<double> qo; // observed q values
    vector<double> Io; // observed I values
    vector<double> Im; // model I values
    vector<double> sigma; // error in Io

    vector<double>& xm; // full model x data
    vector<double>& ym; // full model y data

    /**
     * @brief Calculate chi2 for a given choice of parameters @a params.
     */
    double chi2(const double* params) const {
        double k = params[0];
        double a = params[1];
        double chi = 0;
        for (size_t i = 0; i < qo.size(); i++) {
            double I = k*Im[i] + a;
            chi += pow((Io[i] - I)/sigma[i], 2);
        }
        return chi;
    }

    /**
     * @brief Prepare this class for fitting.
     * @param file measured values to compare the model against.
     * @param x model q values.
     * @param y model I values.
     */
    void setup(string file, vector<double>& x, vector<double>& y) {
        std::tie(qo, Io, sigma) = read(file); // observed values

        // since our model is not analytically derived, we cannot evaluate it at any random point we'd like
        // thus we have to interpolate 
        Im = vector<double>(qo.size()); // model values
        CubicSpline s(x, y);
        for (size_t i = 0; i < qo.size(); ++i) {
            Im[i] = s.spline(qo[i]);
        }
    }

    std::tuple<vector<double>, vector<double>, vector<double>> read(string file) const {
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

            // add the values to our vectors
            q.push_back(_q);
            I.push_back(_I);
            sigma.push_back(_sigma); 
        }
        return std::make_tuple(q, I, sigma);
    }
};