#include <rigidbody/sequencer/ConstraintIteratorElementCallback.h>
#include <rigidbody/sequencer/ConstraintIteratorElement.h>

using namespace rigidbody::sequencer;

ConstraintIteratorElementCallback::ConstraintIteratorElementCallback(ConstraintIteratorElement* caller) : LoopElementCallback(caller->owner), caller(caller) {}