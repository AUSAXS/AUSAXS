// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/OnImprovementElement.h>
#include <rigidbody/detail/Configuration.h>

#include <limits>

using namespace ausaxs::rigidbody::sequencer;

OnImprovementElement::OnImprovementElement(observer_ptr<LoopElement> owner) : LoopElement(owner, 1), best_chi2(std::numeric_limits<double>::max()) {}

OnImprovementElement::~OnImprovementElement() = default;

void OnImprovementElement::run() {
    double new_best_conf = _get_best_conf()->chi2;
    if (new_best_conf < best_chi2) {
        best_chi2 = new_best_conf;
        for (auto& e : elements) {
            e->run();
        }
    }
}