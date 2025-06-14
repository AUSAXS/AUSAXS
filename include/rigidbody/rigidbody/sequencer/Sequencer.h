// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/RigidbodyFwd.h>
#include <rigidbody/sequencer/elements/setup/SetupElement.h>
#include <rigidbody/sequencer/elements/LoopElement.h>
#include <data/DataFwd.h>
#include <utility/observer_ptr.h>

namespace ausaxs::rigidbody::sequencer {
    class Sequencer : public LoopElement {
        public:
            Sequencer();
            Sequencer(const io::ExistingFile& saxs);
            ~Sequencer();

            /**
             * @brief Execute the sequencer.
             * @return The result of the fit.
             */
            std::shared_ptr<fitter::FitResult> execute() override;

            /**
             * @brief Get the Rigidbody object.
             */
            observer_ptr<Rigidbody>& _get_rigidbody();
            observer_ptr<Rigidbody> _get_rigidbody() const override;

            /**
             * @brief Get the molecule object.
             * 
             * This is a convenience method to access the molecule from the Rigidbody.
             */
            observer_ptr<data::Molecule> _get_molecule() const override;

            /**
             * @brief Get the best configuration.
             */
            observer_ptr<detail::BestConf> _get_best_conf() const override;

            /**
             * @brief Get the top Sequencer object.
             */
            observer_ptr<const Sequencer> _get_sequencer() const override;

            observer_ptr<SetupElement> setup();

            /**
             * @brief Perform an optimization step on the rigid body.
             * @return True if a better configuration was found, false otherwise.
             */
            bool _optimize_step() const;

        private:
            SetupElement setup_loop;
            observer_ptr<Rigidbody> rigidbody;
            std::unique_ptr<detail::BestConf> best;
    };
}