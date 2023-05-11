#include <rigidbody/sequencer/ParameterElementCallback.h>
#include <rigidbody/sequencer/ParameterElement.h>

using namespace rigidbody::sequencer;

ParameterElementCallback::ParameterElementCallback(ParameterElement* caller) : LoopElementCallback(caller->owner), caller(caller) {}

ParameterElementCallback::~ParameterElementCallback() = default;

ParameterElement& ParameterElementCallback::amplitude(double amplitude) {
    return caller->amplitude(amplitude);
}