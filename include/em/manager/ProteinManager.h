#pragma once

#include <data/DataFwd.h>
#include <hist/HistFwd.h>

#include <memory>
#include <vector>

namespace em {
    class ImageStackBase;

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
                ProteinManager(const em::ImageStackBase* images);

                /**
                 * @brief Destructor.
                 */
                virtual ~ProteinManager() = default;

                /**
                 * @brief Get the Protein backing this object. 
                 */
                virtual data::Molecule* get_protein() const = 0;

                /**
                 * @brief Get the Protein generated from a given cutoff.
                 */
                virtual data::Molecule* get_protein(double cutoff) = 0;

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
                const em::ImageStackBase* images; 
                std::vector<double> charge_levels;
        };
    }
}