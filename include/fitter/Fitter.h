#pragma once

#include <string>
#include <vector>
#include <map>
#include <iomanip>
#include <iostream>
#include <utility>
#include <memory>

#include <dataset/Multiset.h>
#include <dataset/Dataset.h>
#include <fitter/Fit.h>

class Fitter {
  public:
    virtual ~Fitter() {}

    virtual std::shared_ptr<Fit> fit() = 0;

    /**
     * @brief Make a plot of the fit. 
     * 
     * @return A vector of TGraphs {Interpolated points, Optimal line, Measured points with uncertainties}
     */
    virtual Fit::Plots plot() = 0;

    /**
     * @brief Make a residual plot of the fit.
     * 
     * @return A TGraphErrors with the residuals and their uncertainties. 
     */
    virtual SimpleDataset plot_residuals() = 0;

    virtual std::shared_ptr<Fit> get_fit() const = 0;

    virtual unsigned int dof() const = 0;
};