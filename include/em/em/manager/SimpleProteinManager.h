// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <em/manager/SmartProteinManager.h>

namespace ausaxs::em::managers {
    /**
     * @brief A ProteinManager that regenerates the whole Molecule on every cutoff change.
     *
     * Unlike its base SmartProteinManager, which incrementally replaces only the bodies affected by
     * a new cutoff, this manager rebuilds the entire Molecule from scratch each time. It is simpler
     * and uses a non-partial histogram manager, trading the incremental speed-up for predictability.
     */
    class SimpleProteinManager : public SmartProteinManager {
        public:
            using SmartProteinManager::SmartProteinManager;
            ~SimpleProteinManager() override = default;

        protected:
            /**
             * @brief Rebuild the underlying Molecule from all atoms above the given cutoff.
             */
            void update_protein(double cutoff) override;
    };
}