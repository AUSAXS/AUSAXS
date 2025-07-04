// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hydrate/culling/CullingStrategy.h>

using namespace ausaxs;
using namespace ausaxs::hydrate;

CullingStrategy::CullingStrategy(observer_ptr<data::Molecule> molecule) : molecule(molecule) {}

CullingStrategy::~CullingStrategy() = default;

void CullingStrategy::set_target_count(unsigned int target_count) {this->target_count = target_count;}
