#pragma once

#include <em/detail/EMInternalFwd.h>
#include <data/DataFwd.h>
#include <hist/HistFwd.h>
#include <utility/observer_ptr.h>

#include <memory>
#include <vector>

namespace em {
    namespace managers {
        /**
         * @brief A helper class for the ImageStack. 
         * 
         * This class is responsible for generating and updating histograms
         */
        class ProteinManager {
            public:
                /**
                 * @brief Construct a Manager from an ImageStack.
                 */
                ProteinManager(observer_ptr<const em::ImageStackBase> images);

                /**
                 * @brief Destructor.
                 */
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
                 * @brief Set the charge levels.
                 */
                virtual void set_charge_levels(std::vector<double> levels) noexcept;

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