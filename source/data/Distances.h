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
    Distances(const vector<double>& p_pp, const vector<double>& p_hh, const vector<double>& p_hp, const vector<double>& p_tot)
        : p_pp(p_pp), p_hh(p_hh), p_hp(p_hp), p_tot(p_tot) {setup_bin_dists();}

    void setup_bin_dists();

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

private:
    vector<double> p_pp, p_hh, p_hp, p_tot; // binned distances
    vector<double> d; // the distance corresponding to each bin

    /**
     * @brief Calculate the guinier approximation of the scattering intensity. 
     * @return log10 I(q)
     */
    vector<double> calc_guinier_approx() const;
};