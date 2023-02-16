#pragma once

#include <em/ProteinManager.h>

namespace em {
    //! should probably be moved to the tests where it belongs
    class SimpleProteinManager : public ProteinManager {
        public: 
            using ProteinManager::ProteinManager;

            ~SimpleProteinManager() override = default;

            void update_protein(double cutoff) override;

            void set_charge_levels(std::vector<double> levels) noexcept override;

            void set_charge_levels(Axis levels) noexcept override;
    };
}