// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/HistFwd.h>
#include <hist/detail/BodyTracker.h>
#include <settings/HistogramSettings.h>

#include <memory>

namespace ausaxs::hist {
	/**
	 * @brief A generic histogram manager interface. 
     *        This class is responsible for calculating all distances for a given molecule.
	 */
	class IHistogramManager {
		public:
			virtual ~IHistogramManager() = default;

			/**
			 * @brief Calculate only the total scattering histogram. 
			 */
			virtual std::unique_ptr<DistanceHistogram> calculate() = 0;

			/**
			 * @brief Calculate all contributions to the scattering histogram. 
			 */
			virtual std::unique_ptr<IScatteringModel> calculate_all() = 0;

			/**
			 * @brief The type of this manager, such that the factory can construct another one like it. 
			 */
			virtual settings::hist::HistogramManagerChoice get_type() const = 0;

			/**
			 * @brief The scatterers of the molecule as this manager weighs them, grouped in the layout of the scattering model built by calculate_all(). 
			 *        These are the input of IScatteringModel::amplitude_transform().
			 *        Throws by default, for the managers without an amplitude.
			 */
			virtual AmplitudeScatterers get_amplitude_scatterers() const;
    };
}