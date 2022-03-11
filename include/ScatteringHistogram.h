#pragma once

#include <vector>
#include <string>
#include <utility>
#include "Histogram.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <iostream>

using std::vector, std::string, std::shared_ptr, std::unique_ptr;
using namespace ROOT;

class ScatteringHistogram : Histogram {
  public:
    /**
     * @brief Default constructor.
     */
    ScatteringHistogram() {
        setup();
    }

    /**
     * @brief Move constructor.
     */
    ScatteringHistogram(const ScatteringHistogram&& sh) noexcept : Histogram(sh.p, sh.axis), _p_pp(sh.p_pp), _p_hh(sh.p_hh), _p_hp(sh.p_hp) {
        setup();
    }

    /**
     * @brief Copy constructor. 
     */
    ScatteringHistogram(const ScatteringHistogram& sh) : Histogram(sh.p, sh.axis), _p_pp(sh.p_pp), _p_hh(sh.p_hh), _p_hp(sh.p_hp) {
        setup();
    }

    ScatteringHistogram(const vector<double>& p_pp, const vector<double>& p_hh, const vector<double>& p_hp, const vector<double>& p_tot, const Axis& axis)
        : Histogram(p_tot, axis), _p_pp(p_pp), _p_hh(p_hh), _p_hp(p_hp) {setup();}

    /**
     * @brief Applies the scaling factor @a k to the contribution from the water molecules to this histogram. 
     *        Only affects the total histogram @a p_tot.
     */
    void apply_water_scaling_factor(const double& k);

    /**
     * @brief Removes any scaling applied to the water molecules.
     */
    void reset_water_scaling_factor() {apply_water_scaling_factor(1);}

    /**
     * @brief Prepare a plot of the distances contained in this class.
     * 
     * @return A vector of histograms of the form (atom-atom hist, water-water hist, atom-water hist, total hist))
     */
    vector<shared_ptr<TH1D>> plot_distance() const;

    /**
     * @brief Prepare a plot of the Debye scattering intensities.
     * 
     * @return A histogram of the scattering intensity. 
     */
    unique_ptr<TH1D> plot_debye_scattering() const;

    /**
     * @brief Prepare a plot of the Guinier gyration ratio. 
     * 
     * @return A histogram with logarithmic base-10 y-axis. 
     */
    unique_ptr<TH1D> plot_guinier_approx() const;

    /**
     * @brief Calculate the squared Guinier gyration ratio. 
     */
    double calc_guinier_gyration_ratio_squared() const;

    /**
     * @brief Calculate the intensity based on the Debye scattering equation.
     * 
     * @param r the square of this parameter will be multiplied onto each I(q).
     * 
     * @return I(q)
     */
    vector<double> calc_debye_scattering_intensity() const;

    /**
     * @brief Assign another ScatteringHistogram to this object.
     */
    ScatteringHistogram& operator=(const ScatteringHistogram& h);

    /**
     * @brief Assign another ScatteringHistogram to this object.
     */
    ScatteringHistogram& operator=(ScatteringHistogram&& h);

    const vector<double>& q = _q;
    const vector<double>& p_pp = _p_pp;
    const vector<double>& p_hh = _p_hh;
    const vector<double>& p_hp = _p_hp;
    const vector<double>& p_tot = p;

  private:
    vector<double> _p_pp, _p_hh, _p_hp; // binned distances
    vector<double> _d; // the distance corresponding to each bin
    vector<double> _q; // the q values used as the x-axis

    /**
     * @brief Calculate the guinier approximation of the scattering intensity. 
     * 
     * @return log10 I(q)
     */
    vector<double> calc_guinier_approx() const;

    void setup();
};