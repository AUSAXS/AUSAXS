// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/sequencer/elements/LoopElementCallback.h>
#include <rigidbody/parameters/decay/DecayStrategy.h>

namespace ausaxs::rigidbody::sequencer {
    class ParameterElementCallback : public LoopElementCallback {
        public:
            ParameterElementCallback(ParameterElement* caller);
            virtual ~ParameterElementCallback() override;

            ParameterElement& max_rotation_angle(double radians);

            ParameterElement& max_translation_distance(double distance);

            ParameterElement& decay_strategy(std::unique_ptr<rigidbody::parameter::decay::DecayStrategy> strategy);

        private:
            ParameterElement* caller;
    };
}