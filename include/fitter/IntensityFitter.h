#pragma once

#include "fitter/Fitter.h"
#include "fitter/SimpleLeastSquares.h"
#include "math/CubicSpline.h"
#include "settings.h"

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <tuple>
#include <map>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>

#include "ScatteringHistogram.h"
#include "Exceptions.h"

#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <TGraph.h>
#include <TGraphErrors.h>

using std::string, std::vector, std::shared_ptr, std::unique_ptr;

class IntensityFitter : public Fitter {
  public: 
    /**
     * @brief Prepare a fit of the measured values in @a input to the model described by @a q and @a I.
     * @param input the path to the file containing the measured values. 
     * @param q the model q values.
     * @param I the model I values. 
     */
    // IntensityFitter(string input, vector<double>& q, vector<double>& I) : xm(q), ym(I) {setup(input, q, I);}
    // IntensityFitter(string input, ScatteringHistogram& h) : h(h) {setup(input);}
    IntensityFitter(string input, std::shared_ptr<ScatteringHistogram> h) : h(h), xm(h->q) {setup(input);}
    ~IntensityFitter() override {}

    /**
     * @brief Perform the fit.
     * @return A Fit object containing various information about the fit. Note that the fitted scaling parameter is a = c/M*r_e^2 and b = background
     */
    shared_ptr<Fit> fit() override;

    // shared_ptr<Fit> fit() override {
    //     SimpleLeastSquares fitter(Im_spliced, Io, sigma);
    //     fitted = fitter.fit();
    //     return fitted;
    // }

    vector<shared_ptr<TGraph>> plot() const;

    unique_ptr<TGraphErrors> plot_residuals() const;

  private: 
    shared_ptr<Fit> fitted;
    std::shared_ptr<ScatteringHistogram> h;
    vector<double> qo; // observed q values
    vector<double> Io; // observed I values
    vector<double> sigma; // error in Io

    const vector<double> &xm; // full x and y data from the model

    /**
     * @brief Calculate chi2 for a given choice of parameters @a params.
     */
    double chi2(const double* params) const;

    /**
     * @brief Prepare this class for fitting.
     * @param file measured values to compare the model against.
     */
    void setup(string file);

    /**
     * @brief Splice values from the model to fit the evaluation points defined by the q values of the input file. 
     * @param ym the model y-values corresponding to xm
     */
    vector<double> splice(const vector<double>& ym) const;

    /**
     * @brief Load a data file containing the observed I values. 
     * @param file The file to be read. 
     */
    std::tuple<vector<double>, vector<double>, vector<double>> read(string file) const;
};