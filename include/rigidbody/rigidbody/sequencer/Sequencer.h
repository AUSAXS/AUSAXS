// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/DataFwd.h>
#include <rigidbody/RigidbodyFwd.h>
#include <rigidbody/controller/IController.h>
#include <rigidbody/sequencer/elements/LoopElement.h>
#include <rigidbody/sequencer/elements/setup/SetupElement.h>

#include <memory>

namespace ausaxs::rigidbody::sequencer {
    /**
     * @brief The top-level driver of a rigidbody optimization.
     *
    * The optimization is described declaratively by a sequence script. A Sequencer is itself a
    * LoopElement, and setup() exposes state shared by setup-phase elements. The various `_get_*`
    * accessors give nested elements access to that state. execute() runs the parsed sequence and
    * returns the resulting fit.
     */
    class Sequencer : public LoopElement {
        friend class SetupElement;
        public:
            Sequencer();
            Sequencer(const io::ExistingFile& saxs);
            ~Sequencer() override;

            /**
             * @brief Execute the sequencer.
             * @return The result of the fit.
             */
            std::shared_ptr<fitter::FitResult> execute() override;

            /**
             * @brief Disabled. Sequencer must always be invoked via execute() to ensure proper setup.
             *        Calling run() directly skips rigidbody and controller initialization.
             */
            void run() override;

            /**
             * @brief Expose state shared by setup-phase elements.
             */
            SetupElement& setup();

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