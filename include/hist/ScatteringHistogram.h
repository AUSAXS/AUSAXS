#pragma once

#include <hist/Histogram.h>
#include <hist/DebyeLookupTable.h>
#include <dataset/SimpleDataset.h>

#include <vector>
#include <string>
#include <utility>

namespace hist {
	class ScatteringHistogram : public Histogram {
		public:
			/**
			 * @brief Default constructor.
			 */
			ScatteringHistogram();

			/**
			 * @brief Move constructor.
			 */
			ScatteringHistogram(const ScatteringHistogram&& sh) noexcept;

			/**
			 * @brief Copy constructor. 
			 */
			ScatteringHistogram(const ScatteringHistogram& sh);

			ScatteringHistogram(const std::vector<double>& p_pp, const std::vector<double>& p_hh, const std::vector<double>& p_hp, const std::vector<double>& p_tot, const Axis& axis);

			~ScatteringHistogram() override;

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
			[[nodiscard]] std::vector<Histogram> plot_distance() const;

			/**
			 * @brief Prepare a plot of the Debye scattering intensities.
			 * 
			 * @return A histogram of the scattering intensity. 
			 */
			[[nodiscard]] Histogram plot_debye_scattering() const;

			/**
			 * @brief Prepare a plot of the Guinier gyration ratio. 
			 * 
			 * @return I(q)
			 */
			[[nodiscard]] Histogram plot_guinier_approx() const;

			/**
			 * @brief Calculate the squared Guinier gyration ratio. 
			 */
			[[nodiscard]] double calc_guinier_gyration_ratio_squared() const;

			/**
			 * @brief Calculate the intensity based on the Debye scattering equation.
			 * 
			 * @return I(q)
			 */
			[[nodiscard]] SimpleDataset calc_debye_scattering_intensity() const;

			/**
			 * @brief Calculate the intensity based on the Debye scattering equation for a specific set of scattering vectors.
			 * 		  Note that this recalculates the sin(qd)/qd values for the given q values, and is thus significantly slower than the version without arguments. 
			 * 
			 * @param q The scattering vectors to calculate the intensity for. 
			 * 
			 * @return I(q)
			 */
			[[nodiscard]] SimpleDataset calc_debye_scattering_intensity(const std::vector<double>& q) const;

			ScatteringHistogram& operator=(const ScatteringHistogram& h);
			ScatteringHistogram& operator=(ScatteringHistogram&& h);
			ScatteringHistogram& operator+=(const ScatteringHistogram& rhs);
			ScatteringHistogram& operator-=(const ScatteringHistogram& rhs);
			ScatteringHistogram& operator*=(double rhs);

			bool operator==(const ScatteringHistogram& h) const;

			/**
			 * @brief Extend the view axis to the given maximum value.
			 */
			// void extend_axis(double qmax);

            [[nodiscard]] std::string to_string() const noexcept override;

			Histogram p_pp, p_hh, p_hp; // binned distances
			std::vector<double> d; // The physical distance corresponding to each bin.
			std::vector<double> q; // The q values used as the x-axis.

		private:
			std::unique_ptr<table::DebyeLookupTable> sinqd_table; // Lookup-table for sin(qd)/qd values for the scattering histograms.

			/**
			 * @brief Calculate the guinier approximation of the scattering intensity. 
			 * 
			 * @return log10 I(q)
			 */
			[[nodiscard]] SimpleDataset calc_guinier_approx() const;

			void setup();
	};

	ScatteringHistogram operator+(const ScatteringHistogram& lhs, const ScatteringHistogram& rhs);
	ScatteringHistogram operator-(const ScatteringHistogram& lhs, const ScatteringHistogram& rhs);
	ScatteringHistogram operator*(const ScatteringHistogram& lhs, double rhs);
}