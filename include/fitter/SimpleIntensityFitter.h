#pragma once

#include "fitter/Fitter.h"
#include "math/SimpleLeastSquares.h"
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

class SimpleIntensityFitter : public Fitter {
  public: 
    /**
     * @brief Constructor.
     * 
     * Prepare a fit of the measured values in @a input to a model to be defined later. 
     * 
     * @param input The path to the file containing the measured values. 
     */
    SimpleIntensityFitter(string input) {setup(input);}

    /**
     * @brief Constructor.
     * 
     * Prepare a fit of the measured values in @a input to the model described by @a h.
     * 
     * @param input The path to the file containing the measured values. 
     * @param h The ScatteringHistogram to fit. 
     */
    SimpleIntensityFitter(string input, const ScatteringHistogram& h) : h(h) {setup(input);}

    /**
     * @brief Constructor.
     * 
     * Prepare a fit of the measured values in @a input to the model described by @a h.
     * 
     * @param input the path to the file containing the measured values. 
     * @param h The ScatteringHistogram to fit. 
     */
    SimpleIntensityFitter(string input, ScatteringHistogram&& h) : h(std::move(h)) {setup(input);}

    /**
     * @brief Constructor. 
     * 
     * Prepare a fit to the data set defined by @a q, @a I, and @a sigma. 
     * 
     * @param q The scattering vectors.  
     * @param I The corresponding intensities.  
     * @param sigma The uncertainties in each measurement. 
     */
    SimpleIntensityFitter(vector<double>& q, vector<double>& I, vector<double>& sigma) : qo(q), Io(I), sigma(sigma) {}

    /**
     * @brief Constructor. 
     * 
     * Prepare a fit of the histogram @a h to the data set defined by @a q, @a I, and @a sigma. 
     * 
     * @param q The scattering vectors.  
     * @param I The corresponding intensities.  
     * @param sigma The uncertainties in each measurement. 
     * @param h The histogram to be fitted. 
     */
    SimpleIntensityFitter(vector<double>& q, vector<double>& I, vector<double>& sigma, const ScatteringHistogram& h) : qo(q), Io(I), sigma(sigma), h(h) {}

    /**
     * @brief Constructor.
     * 
     * Prepare a fit of the first histogram to the second. A series of data points is extracted from @a h2 and used to fit @a h1.
     * 
     * @param data The data histogram. 
     * @param model The model histogram. 
     * @param limits The limits on the generated data points. 
     */
    SimpleIntensityFitter(const ScatteringHistogram& data, const ScatteringHistogram& model, const Limit& limits = Limit(setting::fit::q_low, setting::fit::q_high));

    /**
     * @brief Constructor.
     * 
     * Prepare a fit to the histogram. A series of data points is extracted from it and used as the data points of the model. 
     * 
     * @param model The model histogram. 
     * @param limits The limits on the generated data points. 
     */
    SimpleIntensityFitter(const ScatteringHistogram& model, const Limit& limits = Limit(setting::fit::q_low, setting::fit::q_high));

    /**
     * @brief Destructor.
     */
    ~SimpleIntensityFitter() override = default;

    /**
     * @brief Perform the fit.
     * 
     * @return A Fit object containing various information about the fit. Note that the fitted scaling parameter is a = c/M*r_e^2 and b = background
     */
    virtual shared_ptr<Fit> fit() override;

    /**
     * @brief Make a plot of the fit. 
     * 
     * @return A vector of TGraphs {Interpolated points, Optimal line, Measured points with uncertainties}
     */
    virtual vector<shared_ptr<TGraph>> plot();

    /**
     * @brief Make a residual plot of the fit.
     * 
     * @return A TGraphErrors with the residuals and their uncertainties. 
     */
    virtual unique_ptr<TGraphErrors> plot_residuals();

    /**
     * @brief Change the scattering histogram used for the fit. 
     */
    void set_scattering_hist(const ScatteringHistogram& h);

    /**
     * @brief Change the scattering histogram used for the fit. 
     */
    void set_scattering_hist(ScatteringHistogram&& h);

    /**
     * @brief Get the number of degrees of freedom. 
     */
    unsigned int degrees_of_freedom() const;

    /**
     * @brief Get the number of degrees of freedom. 
     */
    unsigned int dof() const;

  protected: 
    shared_ptr<Fit> fitted;
    vector<double> qo; // observed q values
    vector<double> Io; // observed I values
    vector<double> sigma; // error in Io
    ScatteringHistogram h;

    /**
     * @brief Calculate chi2 for a given choice of parameters @a params.
     */
    virtual double chi2(const double* params);

    /**
     * @brief Prepare this class for fitting.
     * 
     * @param file measured values to compare the model against.
     */
    void setup(string file);

    /**
     * @brief Splice values from the model to fit the evaluation points defined by the q values of the input file. 
     * 
     * @param ym the model y-values corresponding to xm
     */
    vector<double> splice(const vector<double>& ym) const;

    /**
     * @brief Load a data file containing the observed I values. 
     * 
     * @param file The file to be read. 
     */
    std::tuple<vector<double>, vector<double>, vector<double>> read(string file) const;

    /**
     * @brief Initialize this class based on a model histogram. 
     */
    void model_setup(const ScatteringHistogram& model, const Limit& limits);
};