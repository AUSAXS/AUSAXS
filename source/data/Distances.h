#pragma once

// forwards declaration
class Protein;

#include <vector>
#include <string>
#include <utility>
#include "TH1D.h"
#include "TCanvas.h"
#include "Protein.h"

using std::vector, std::string, std::shared_ptr, std::unique_ptr;
using namespace ROOT;

class Distances {
public:
    Distances(Protein* protein, vector<double>& pp, vector<double>& hh, vector<double>& hp, vector<double>& wpp, vector<double>& whh, vector<double>& whp) 
        : protein(protein), d_pp(pp), d_hh(hh), d_hp(hp), w_pp(wpp), w_hh(whh), w_hp(whp), axes(setting::axes::scattering_intensity_plot_binned) {}

    /**
     * @brief Set the plot axes. 
     * @param axes the axes to use in the format {bins, xmin, xmax}.
     */
    void set_axes(const vector<int> axes);

    /**
     * @brief Prepare a plot of the distances contained in this class.
     * @return A vector of histograms of the form (atom-atom hist, water-water hist, atom-water hist, total hist))
     */
    vector<shared_ptr<TH1D>> plot_distance() const;

    /**
     * @brief Prepare a plot of the Debye scattering intensities.
     * @return A histogram of the scattering intensity. 
     */
    unique_ptr<TH1D> plot_debye_scattering() const;

    /**
     * @brief Prepare a plot of the Guinier gyration ratio. 
     * @return A histogram with logarithmic base-10 y-axis. 
     */
    unique_ptr<TH1D> plot_guinier_approx() const;

    unique_ptr<TH1D> fit_debye_plot() const;

    /**
     * @brief Calculate the squared Guinier gyration ratio. 
     */
    double calc_guinier_gyration_ratio_squared() const;

    /**
     * @brief Calculate the intensity based on the Debye scattering equation
     * @param r the square of this parameter will be multiplied onto each I(q).
     * @return I(q)
     */
    vector<double> calc_debye_scattering_intensity() const;

    /**
     * @brief Returns the vector representing the x-axis of q-values. 
     */
    vector<double> get_xaxis() const;

    const vector<double> d_pp, d_hh, d_hp; // raw distances
private:
    const vector<double> w_pp, w_hh, w_hp; // weights on each distance
    vector<int> axes; // bin style in the form {bins, xmin, xmax}
    vector<double> binned_pp, binned_hh, binned_hp, binned_tot; // binned distances
    const Protein* protein; // the protein which created this object

    /**
     * @brief Bin the raw distance data. 
     * @param axes a vector of the form {bins, xmin, xmax} to use for the binning
     */
    void bin(const vector<int>& axes);

    /**
     * @brief Calculate the guinier approximation of the scattering intensity. 
     * @return log10 I(q)
     */
    vector<double> calc_guinier_approx() const;
};