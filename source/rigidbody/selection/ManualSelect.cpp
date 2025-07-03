// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/selection/ManualSelect.h>
#include <utility/Exceptions.h>

using namespace ausaxs::rigidbody::selection;

ManualSelect::ManualSelect(observer_ptr<const RigidBody> rigidbody) : BodySelectStrategy(rigidbody) {}

ManualSelect::~ManualSelect() = default;

std::pair<unsigned int, int> ManualSelect::next() {
    throw except::not_implemented("ManualSelect::next: Not implemented.");
}
