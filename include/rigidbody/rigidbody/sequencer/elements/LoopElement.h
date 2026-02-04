// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/RigidbodyFwd.h>
#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/sequencer/elements/GenericElement.h>
#include <utility/observer_ptr.h>
#include <fitter/FitterFwd.h>
#include <data/DataFwd.h>
#include <io/IOFwd.h>

#include <memory>
#include <vector>

namespace ausaxs::rigidbody::sequencer {
    /**
     * @brief A loop element is a sequence element that repeats whatever is inside it a number of times.
     */
    class LoopElement : public GenericElement {
        friend class OptimizeStepElement;
        public:
            LoopElement(observer_ptr<LoopElement> owner, unsigned int repeats);
            virtual ~LoopElement();

            virtual std::shared_ptr<fitter::FitResult> execute();

            /**
             * @brief Create a nested loop.
             */
            LoopElement& loop(unsigned int repeats);

            /**
             * @brief Set the parameter strategy.
             */
            ParameterElement& parameter_strategy(std::unique_ptr<rigidbody::parameter::ParameterGenerationStrategy> strategy);

            /**
             * @brief Set the body selection strategy.
             */
            BodySelectElement& body_select_strategy(std::unique_ptr<rigidbody::selection::BodySelectStrategy> strategy);

            /**
             * @brief Set the transformation strategy.
             */
            TransformElement& transform_strategy(std::unique_ptr<rigidbody::transform::TransformStrategy> strategy);

            /**
             * @brief Perform a single optimization step.
             */
            OptimizeStepElement& optimize();

            /**
             * @brief End the current loop.
             */
            virtual LoopElement& end();

            /**
             * @brief Save the current state of the system.
             *
             * @param path The path to save the state to. The extension of the file will determine the format.
             */
            LoopElement& save(const io::File& path);

            /**
             * @brief Perform the subroutines for every n iterations of this loop.
             */
            EveryNStepElement& every(unsigned int n);

            /**
             * @brief Run an iteration of this loop. 
             */
            void run() override;

            virtual observer_ptr<Rigidbody> _get_rigidbody() const;

            virtual observer_ptr<data::Molecule> _get_molecule() const;

            virtual observer_ptr<detail::MoleculeTransformParametersAbsolute> _get_best_conf() const;

            virtual observer_ptr<const Sequencer> _get_sequencer() const;
            virtual observer_ptr<Sequencer> _get_sequencer();

            std::vector<std::unique_ptr<GenericElement>>& _get_elements();

            observer_ptr<LoopElement> _get_owner() const;

        protected: 
            unsigned int iterations = 1;
            std::vector<std::unique_ptr<GenericElement>> elements;

        private:
            observer_ptr<LoopElement> owner;
            inline static unsigned int total_loop_count = 0;
            inline static unsigned int global_counter = 0;
    };
}