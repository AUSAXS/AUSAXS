// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <io/IOFwd.h>
#include <rigidbody/RigidbodyFwd.h>
#include <rigidbody/sequencer/SequencerFwd.h>
#include <settings/RigidBodySettings.h>

namespace ausaxs::rigidbody::sequencer {
    /**
        * @brief A callback class for the LoopElement class.
        * 
        * This class is used to provide a callback interface for the LoopElement class, giving other elements access
        * to the basic loop methods.
        */
    class LoopElementCallback {
        public: 
            /**
             * @brief Create a callback interface.
             */
            LoopElementCallback(LoopElement* caller);

            virtual ~LoopElementCallback();

            /**
             * @brief Create a nested loop.
             */
            LoopElement& loop(int repeats) const;

            /**
             * @brief Set the parameter strategy.
             */
            ParameterElement& parameter_strategy(std::unique_ptr<rigidbody::parameter::ParameterGenerationStrategy> strategy) const;

            /**
             * @brief Set the body selection strategy.
             */
            BodySelectElement& body_select_strategy(std::unique_ptr<rigidbody::selection::BodySelectStrategy> strategy) const;

            /**
             * @brief Set the transformation strategy.
             */
            TransformElement& transform_strategy(std::unique_ptr<rigidbody::transform::TransformStrategy> strategy) const;

            /**
             * @brief End the current loop.
             */
            virtual LoopElement& end();

            /**
             * @brief Perform the subroutines for every n iterations of this loop.
             */
            EveryNStepElement& every(int n) const;

            /**
             * @brief Save the current state of the system.
             *
             * @param path The path to save the state to. The extension of the file will determine the format.
             */
            LoopElement& save(const io::File& path) const;

            LoopElement* owner;
    };
}