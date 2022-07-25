#pragma once

#include <vector>
#include <string>
#include <utility>

#include <TH1D.h>

#include <histogram/Histogram.h>
#include <histogram/DebyeLookupTable.h>
#include <utility/SimpleDataset.h>

namespace hist {
	class ScatteringHistogram : public Histogram {
		public:
			/**
			 * @brief Default constructor.
			 */
			ScatteringHistogram() {}

			/**
			 * @brief Move constructor.
			 */
			ScatteringHistogram(const ScatteringHistogram&& sh) noexcept : Histogram(sh.p, sh.axis), p_pp(sh.p_pp), p_hh(sh.p_hh), p_hp(sh.p_hp) {
				setup();
			}

			/**
			 * @brief Copy constructor. 
			 */
			ScatteringHistogram(const ScatteringHistogram& sh) : Histogram(sh.p, sh.axis), p_pp(sh.p_pp), p_hh(sh.p_hh), p_hp(sh.p_hp) {
				setup();
			}

			ScatteringHistogram(const std::vector<double>& p_pp, const std::vector<double>& p_hh, const std::vector<double>& p_hp, const std::vector<double>& p_tot, const Axis& axis)
				: Histogram(p_tot, axis), p_pp(p_pp), p_hh(p_hh), p_hp(p_hp) {setup();}

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
			std::vector<std::shared_ptr<TH1D>> plot_distance() const;

			/**
			 * @brief Prepare a plot of the Debye scattering intensities.
			 * 
			 * @return A histogram of the scattering intensity. 
			 */
			std::unique_ptr<TH1D> plot_debye_scattering() const;

			/**
			 * @brief Prepare a plot of the Guinier gyration ratio. 
			 * 
			 * @return A histogram with logarithmic base-10 y-axis. 
			 */
			std::unique_ptr<TH1D> plot_guinier_approx() const;

			/**
			 * @brief Calculate the squared Guinier gyration ratio. 
			 */
			double calc_guinier_gyration_ratio_squared() const;

			/**
			 * @brief Calculate the intensity based on the Debye scattering equation.
			 * 
			 * @return I(q)
			 */
			SimpleDataset calc_debye_scattering_intensity() const;

			/**
			 * @brief Calculate the intensity based on the Debye scattering equation for a specific set of scattering vectors.
			 * 
			 * @param q The scattering vectors to calculate the intensity for. 
			 * 
			 * @return I(q)
			 */
			SimpleDataset calc_debye_scattering_intensity(const std::vector<double>& q) const;

			/**
			 * @brief Assign another ScatteringHistogram to this object.
			 */
			ScatteringHistogram& operator=(const ScatteringHistogram& h);

			/**
			 * @brief Assign another ScatteringHistogram to this object.
			 */
			ScatteringHistogram& operator=(ScatteringHistogram&& h);

			std::vector<double> p_pp, p_hh, p_hp; // binned distances
			std::vector<double> d; // The distance corresponding to each bin.
			std::vector<double> q; // The q values used as the x-axis.

		private:
			table::DebyeLookupTable sinqd_table; // Lookup-table for sin(qd)/qd values for the scattering histograms.

			/**
			 * @brief Calculate the guinier approximation of the scattering intensity. 
			 * 
			 * @return log10 I(q)
			 */
			SimpleDataset calc_guinier_approx() const;

			void setup();
	};
}