// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/elements/ConvertSymmetryElement.h>

#include <cassert>

using namespace ausaxs::rigidbody::sequencer;

ConvertSymmetryElement::ConvertSymmetryElement(observer_ptr<LoopElement> owner, int) : LoopElementCallback(owner), GenericElement() {
    assert(false && "ConvertSymmetryElement::run: Not implemented.");
}

ConvertSymmetryElement::~ConvertSymmetryElement() = default;

void ConvertSymmetryElement::run() {
    assert(false && "ConvertSymmetryElement::run: Not implemented.");
}