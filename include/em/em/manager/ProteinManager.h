// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <em/detail/EMInternalFwd.h>
#include <data/DataFwd.h>
#include <hist/HistFwd.h>
#include <utility/observer_ptr.h>

#include <memory>
#include <vector>

namespace ausaxs::em {
    namespace managers {
        /**
         * @brief A helper class for the ImageStack. 
         * 
         * This class is responsible for generating and updating histograms
         */
        class ProteinManager {
            public:
                ProteinManager(observer_ptr<const em::ImageStackBase> images);
                virtual ~ProteinManager();

                /**
                 * @brief Get the Protein backing this object. 
                 */
                virtual observer_ptr<const data::Molecule> get_protein() const = 0;

                /**
                 * @brief Get the Protein generated from a given cutoff.
                 */
                virtual observer_ptr<data::Molecule> get_protein(double cutoff) = 0;

                /**
                 * @brief Get the histogram for a given cutoff.
                 */
                virtual std::unique_ptr<hist::ICompositeDistanceHistogram> get_histogram(double cutoff) = 0;

                /**
                 * @brief Get the excluded volume mass of the managed molecule.
                 * 		  This is the excluded volume of the molecule times the average protein mass density. 
                 * 
                 * @return The excluded volume mass in Da.
                 */
                double get_excluded_volume_mass() const;

                /**
                 * @brief Calculate the volume of the managed molecule based on the number of grid bins it spans.
                 * 
                 * @return The volume in Å^3.
                 */
                double get_volume_grid() const;

                /**
                 * @brief Set the charge levels.
                 */
                virtual void set_charge_levels(const std::vector<double>& levels) noexcept;

                /**
                 * @brief Get the charge levels.
                 */
                std::vector<double> get_charge_levels() const noexcept;

            protected:
                observer_ptr<const em::ImageStackBase> images; 
                std::vector<double> charge_levels;
        };
    }
}