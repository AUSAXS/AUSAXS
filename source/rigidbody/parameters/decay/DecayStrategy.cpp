// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/parameters/decay/DecayStrategy.h>
#include <utility/Console.h>

#include <string>

using namespace ausaxs::rigidbody::parameter::decay;

unsigned int DecayStrategy::next_draw() {
    if (iterations != 0 && draws == iterations && !overdrawn_warning_issued) {
        overdrawn_warning_issued = true;
        console::print_warning(
            "DecayStrategy: a parameter generator declared with " + std::to_string(iterations) + " iterations is being "
            "used for more steps than that. Its amplitudes no longer follow the declared decay. Set \"iterations\" to "
            "the number of optimisation steps the generator is actually used for."
        );
    }
    return draws++;
}
