#pragma once

#include <em/ProteinManager.h>

namespace em {
    class SimpleProteinManager : public ProteinManager {
        public: 
            using ProteinManager::ProteinManager;

            ~SimpleProteinManager() override = default;

            /**
             * @brief Update the Protein to reflect a new cutoff value.
             */
            void update_protein(double cutoff) override;

            void set_charge_levels(std::vector<double> levels) noexcept override;

            void set_charge_levels(Axis levels) noexcept override;
    };
}