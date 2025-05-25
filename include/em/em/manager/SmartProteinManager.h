#pragma once

#include <data/DataFwd.h>
#include <hist/HistFwd.h>
#include <em/manager/ProteinManager.h>
#include <em/detail/EMAtom.h>
#include <utility/observer_ptr.h>

#include <vector>

namespace ausaxs::em::managers {
    /**
     * @brief A helper class for the ImageStack. 
     *
     * This class stores a single instance of a Molecule, which is generated from the EM images.
     * Various methods are provided to update the state of this single instance to reflect new cutoff values. 
     */
    class SmartProteinManager : public ProteinManager {
        public:
            using ProteinManager::ProteinManager;
            virtual ~SmartProteinManager() override;

            /**
             * @brief Evaluate and get the histogram for the given cutoff. 
             *        This will update the underlying Molecule to reflect the new cutoff value.
             */
            std::unique_ptr<hist::ICompositeDistanceHistogram> get_histogram(double cutoff) override;

            /**
             * @brief Get the current state of the underlying Molecule. 
             */
            observer_ptr<const data::Molecule> get_protein() const override;

            /**
             * @brief Refresh the underlying Molecule with a new cutoff, and get a pointer to it. 
             */
            observer_ptr<data::Molecule> get_protein(double cutoff) override;

            /**
             * @brief Set the charge levels.
             *        This will invalidate the underlying Molecule, meaning it will be regenerated on the next call to get_protein.
             */
            virtual void set_charge_levels(const std::vector<double>& levels) noexcept override;

        protected:
            std::unique_ptr<data::Molecule> protein;

            /**
             * @brief Generate a list of all atoms with charge densities larger than the cutoff. 
             */
            std::vector<data::EMAtom> generate_atoms(double cutoff) const;

            /**
             * @brief Update the underlying Molecule to reflect a new cutoff value.
             */
            virtual void update_protein(double cutoff);

            /**
             * @brief Enable or disable the histogram manager initialization for generated molecules.
             * 
             * This is used to prevent expensive initialization of the histogram managers whenever a new 
             * molecule is generated, since they will be cannibalized by the member variable anyway. 
             */
            void toggle_histogram_manager_init(bool state);

        private:
            double previous_cutoff = 0;

            /**
             * @brief Generate a new Molecule with the given cutoff. 
             */
            std::unique_ptr<data::Molecule> generate_new_protein(double cutoff) const;
    };
}