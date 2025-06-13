// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/elements/GenericElement.h>
#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/RigidbodyFwd.h>
#include <utility/observer_ptr.h>

namespace ausaxs::rigidbody::sequencer {
    /**
     * @brief Use an existing RigidBody for the optimization. 
     */
    class LoadExistingElement : public GenericElement {
        public:
            /**
             * @brief Use the given RigidBody for the optimization.
             *        The RigidBody must be fully initialized, and must have a lifetime that exceeds that of the Sequencer.
             */
            LoadExistingElement(observer_ptr<Sequencer> owner, observer_ptr<RigidBody> rigidbody);
            virtual ~LoadExistingElement() = default;

            void run() override;
        
        private:
            observer_ptr<Sequencer> owner;
            observer_ptr<RigidBody> rigidbody;
    };
}