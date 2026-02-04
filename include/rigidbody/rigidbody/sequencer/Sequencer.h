// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/RigidbodyFwd.h>
#include <rigidbody/sequencer/elements/setup/SetupElement.h>
#include <rigidbody/sequencer/elements/LoopElement.h>
#include <rigidbody/controller/IController.h>
#include <data/DataFwd.h>

#include <memory>

namespace ausaxs::rigidbody::sequencer {
    class Sequencer : public LoopElement {
        friend class SetupElement;
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
             * @brief Expose the setup element for configuration.
             */
            SetupElement& setup();
            
            LoopElement& end() override;

            /**
             * @brief Get the Rigidbody object.
             */
            observer_ptr<Rigidbody> _get_rigidbody() const override;
            void _set_rigidbody(observer_ptr<Rigidbody> rigidbody);

            /**
             * @brief Get the molecule object.
             * 
             * This is a convenience method to access the molecule from the Rigidbody.
             */
            observer_ptr<data::Molecule> _get_molecule() const override;

            /**
             * @brief Get the top Sequencer object.
             */
            observer_ptr<const Sequencer> _get_sequencer() const override;
            observer_ptr<Sequencer> _get_sequencer() override;

            /**
             * @brief Get the current rigidbody controller.
             */
            observer_ptr<controller::IController> _get_controller() const;

            /**
             * @brief Get the best configuration found so far.
             */
            observer_ptr<rigidbody::detail::MoleculeTransformParametersAbsolute> _get_best_conf() const override;
            
        private:
            SetupElement setup_loop;
            observer_ptr<Rigidbody> rigidbody;
    };
}