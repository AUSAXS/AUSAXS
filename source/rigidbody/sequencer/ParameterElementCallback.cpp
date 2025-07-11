// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/ParameterElementCallback.h>
#include <rigidbody/sequencer/ParameterElement.h>

using namespace ausaxs::rigidbody::sequencer;

ParameterElementCallback::ParameterElementCallback(ParameterElement* caller) : LoopElementCallback(caller->owner), caller(caller) {}

ParameterElementCallback::~ParameterElementCallback() = default;

ParameterElement& ParameterElementCallback::max_rotation_angle(double radians) {
    return caller->max_rotation_angle(radians);
}

ParameterElement& ParameterElementCallback::max_translation_distance(double distance) {
    return caller->max_translation_distance(distance);
}

ParameterElement& ParameterElementCallback::decay_strategy(std::unique_ptr<rigidbody::parameter::decay::DecayStrategy> strategy) {
    return caller->decay_strategy(std::move(strategy));
}