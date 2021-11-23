#pragma once

#include <vector>
#include <string>
#include <utility>
#include "TH1D.h"
#include "TCanvas.h"

using std::vector, std::string, std::shared_ptr, std::unique_ptr;
using namespace ROOT;

class Distances {
public:
    Distances(vector<double>& pp, vector<double>& hh, vector<double>& hp, vector<double>& wpp, vector<double>& whh, vector<double>& whp) 
        : d_pp(pp), d_hh(hh), d_hp(hp), w_pp(wpp), w_hh(whh), w_hp(whp) {}

    /**
     * @brief Prepare a plot of the distances contained in this class.
     * @param axes the axes to use in the format {bins, xmin, xmax}
     * @return A pair (canvas, vector(atom-atom hist, water-water hist, atom-water hist, total hist))
     */
    std::pair<unique_ptr<TCanvas>, vector<shared_ptr<TH1D>>> plot_distance(const vector<int>& axes);

    /**
     * @brief Prepare a plot of the Debye scattering intensities.
     * @return A pair (canvas, hist)
     */
    std::pair<unique_ptr<TCanvas>, unique_ptr<TH1D>> plot_debye_scattering();

    /**
     * @brief Prepare a plot of the Guinier gyration ratio. 
     * @return A pair (canvas, hist)
     */
    std::pair<unique_ptr<TCanvas>, unique_ptr<TH1D>> plot_Guinier_gyration();

    const vector<double> d_pp, d_hh, d_hp; // raw distances
private:
    const vector<double> w_pp, w_hh, w_hp; // weights on each distance
    vector<int> axes; // bin style in the form {bins, xmin, xmax}
    vector<double> binned_pp, binned_hh, binned_hp, binned_tot; // binned distances

    /**
     * @brief Bin the raw distance data. 
     * @param axes a vector of the form {bins, xmin, xmax} to use for the binning
     */
    void bin(const vector<int>& axes);

    /**
     * @brief Calculate the intensity based on the Debye scattering equation
     */
    vector<double> calc_debye_scattering_intensity() const;

    /**
     * @brief Calculate the Guinier gyration ratio
     */
    void calc_gyration_ratio() const;
};