// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/ConstraintIteratorElementCallback.h>
#include <rigidbody/sequencer/elements/ConstraintIteratorElement.h>

using namespace ausaxs::rigidbody::sequencer;

ConstraintIteratorElementCallback::ConstraintIteratorElementCallback(ConstraintIteratorElement* caller) : LoopElementCallback(caller->owner), caller(caller) {}