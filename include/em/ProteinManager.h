#pragma once

#include <vector>

#include <em/detail/ImageStackBase.h>
#include <em/CullingStrategy.h>
#include <data/Protein.h>
#include <hist/ScatteringHistogram.h>

namespace em {
    /**
     * @brief A helper class for the ImageStack. 
     * 
     * This class generates and updates histograms in a smart way. This is done by splitting the 
     * generated atoms into multiple bodies, and then utilizing the smart histogram manager from the Protein. 
     */
    class ProteinManager {
        public:
            /**
             * @brief Construct a Manager from an ImageStack.
             */
            ProteinManager(const ImageStackBase& images);

            /**
             * @brief Destructor.
             */
            virtual ~ProteinManager() = default;

            /**
             * @brief Get the histogram for a given cutoff.
             */
            hist::ScatteringHistogram get_histogram(double cutoff);

            /**
             * @brief Get the Protein backing this object. 
             */
            std::shared_ptr<Protein> get_protein() const;

            /**
             * @brief Get the Protein generated from a given cutoff.
             */
            std::shared_ptr<Protein> get_protein(double cutoff);

            /**
             * @brief Set the charge levels to default values.
             */
            virtual void set_charge_levels() noexcept;

            /**
             * @brief Set the charge levels.
             */
            virtual void set_charge_levels(std::vector<double> levels) noexcept;

            /**
             * @brief Set the charge levels.
             */
            virtual void set_charge_levels(Axis levels) noexcept;

            std::vector<double> get_charge_levels() const noexcept;

        protected:
            const ImageStackBase& images; 
            std::shared_ptr<Protein> protein;

            /**
             * @brief Generate the atmos for a given cutoff.
             */
            std::vector<Atom> generate_atoms(double cutoff) const;

            /**
             * @brief Update the Protein to reflect a new cutoff value.
             */
            virtual void update_protein(double cutoff);

        private:
            std::vector<double> charge_levels;
            double previous_cutoff = 0;

            /**
             * @brief Generate a new Protein for a given cutoff. 
             */
            std::unique_ptr<Protein>generate_protein(double cutoff) const;
    };
}