#pragma once

#include <em/manager/ProteinManager.h>
#include <data/DataFwd.h>
#include <hist/HistFwd.h>
#include <utility/observer_ptr.h>

#include <vector>

namespace em::managers {
    /**
     * @brief A helper class for the ImageStack. 
     * 
     * This class generates and updates histograms in a smart way. This is done by splitting the 
     * generated atoms into multiple bodies, and then utilizing the smart histogram manager from the Protein. 
     */
    class SmartProteinManager : public ProteinManager {
        public:
            using ProteinManager::ProteinManager;

            /**
             * @brief Destructor.
             */
            virtual ~SmartProteinManager() = default;

            /**
             * @brief Get the histogram for a given cutoff.
             */
            std::unique_ptr<hist::ICompositeDistanceHistogram> get_histogram(double cutoff) override;

            /**
             * @brief Get the Protein backing this object. 
             */
            observer_ptr<const data::Molecule> get_protein() const override;

            /**
             * @brief Get the Protein generated from a given cutoff.
             */
            observer_ptr<data::Molecule> get_protein(double cutoff) override;

            /**
             * @brief Set the charge levels.
             */
            virtual void set_charge_levels(const std::vector<double>& levels) noexcept;

        protected:
            std::unique_ptr<data::Molecule> protein;

            /**
             * @brief Generate the atmos for a given cutoff.
             */
            std::vector<data::record::Atom> generate_atoms(double cutoff) const;

            /**
             * @brief Update the Protein to reflect a new cutoff value.
             */
            virtual void update_protein(double cutoff);

        private:
            double previous_cutoff = 0;

            /**
             * @brief Generate a new Protein for a given cutoff. 
             */
            std::unique_ptr<data::Molecule> generate_protein(double cutoff) const;
    };
}