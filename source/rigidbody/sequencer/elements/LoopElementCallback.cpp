// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/LoopElementCallback.h>

#include <rigidbody/sequencer/elements/BodySelectElement.h>
#include <rigidbody/sequencer/elements/LoopElement.h>
#include <rigidbody/sequencer/elements/ParameterElement.h>
#include <rigidbody/sequencer/elements/TransformElement.h>

using namespace ausaxs::rigidbody::sequencer;

LoopElementCallback::LoopElementCallback(LoopElement* caller) : owner(caller) {}

LoopElementCallback::~LoopElementCallback() = default;

LoopElement& LoopElementCallback::loop(int repeats) const {
    return owner->loop(repeats);
}

ParameterElement& LoopElementCallback::parameter_strategy(std::unique_ptr<rigidbody::parameter::ParameterGenerationStrategy> strategy) const {
    return owner->parameter_strategy(std::move(strategy));
}

BodySelectElement& LoopElementCallback::LoopElementCallback::body_select_strategy(std::unique_ptr<rigidbody::selection::BodySelectStrategy> strategy) const {
    return owner->body_select_strategy(std::move(strategy));
}

TransformElement& LoopElementCallback::transform_strategy(std::unique_ptr<rigidbody::transform::TransformStrategy> strategy) const {
    return owner->transform_strategy(std::move(strategy));
}

LoopElement& LoopElementCallback::end() {
    return owner->end();
}

LoopElement& LoopElementCallback::save(const io::File& path) const {
    return owner->save(path);
}

EveryNStepElement& LoopElementCallback::every(int n) const {
    return owner->every(n);
}