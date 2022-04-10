#pragma once

#include <string>
#include <vector>
#include <map>
#include <iomanip>
#include <iostream>
#include <utility>
#include <memory>

#include <Math/Minimizer.h>
#include <TGraph.h>
#include <TGraphErrors.h>

class Fit;

class Fitter {
  public:
    virtual ~Fitter() {}

    virtual std::shared_ptr<Fit> fit() = 0;

    /**
     * @brief Make a plot of the fit. 
     * 
     * @return A vector of TGraphs {Interpolated points, Optimal line, Measured points with uncertainties}
     */
    virtual std::vector<std::shared_ptr<TGraph>> plot() = 0;

    /**
     * @brief Make a residual plot of the fit.
     * 
     * @return A TGraphErrors with the residuals and their uncertainties. 
     */
    virtual std::unique_ptr<TGraphErrors> plot_residuals() = 0;

    virtual unsigned int dof() const = 0;
};