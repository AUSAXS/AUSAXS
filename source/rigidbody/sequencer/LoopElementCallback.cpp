/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/sequencer/LoopElementCallback.h>
#include <rigidbody/sequencer/LoopElement.h>
#include <rigidbody/sequencer/BodySelectElement.h>
#include <rigidbody/sequencer/TransformElement.h>
#include <rigidbody/sequencer/ParameterElement.h>

using namespace rigidbody::sequencer;

LoopElementCallback::LoopElementCallback(LoopElement* caller) : owner(caller) {}

LoopElementCallback::~LoopElementCallback() = default;

LoopElement& LoopElementCallback::loop(unsigned int repeats) {
    return owner->loop(repeats);
}

ParameterElement& LoopElementCallback::parameter_strategy(std::unique_ptr<rigidbody::parameter::ParameterGenerationStrategy> strategy) {
    return owner->parameter_strategy(std::move(strategy));
}

BodySelectElement& LoopElementCallback::LoopElementCallback::body_select_strategy(std::unique_ptr<rigidbody::selection::BodySelectStrategy> strategy) {
    return owner->body_select_strategy(std::move(strategy));
}

TransformElement& LoopElementCallback::transform_strategy(std::unique_ptr<rigidbody::transform::TransformStrategy> strategy) {
    return owner->transform_strategy(std::move(strategy));
}

LoopElement& LoopElementCallback::end() {
    return owner->end();
}

LoopElement& LoopElementCallback::save(const io::File& path) {
    return owner->save(path);
}

EveryNStepElement& LoopElementCallback::every(unsigned int n) {
    return owner->every(n);
}