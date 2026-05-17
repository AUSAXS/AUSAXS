// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <data/symmetry/ISymmetry.h>

using namespace ausaxs::symmetry;

std::vector<CopyPair> ISymmetry::internal_pair_schedule() const {
    // Cyclic-chain reuse. With R = repetitions(), the bodies are {original, copy_1, ..., copy_R}.
    // Copy k sits one fixed generator step from copy k-1, so the distance between two bodies
    // depends only on their index separation s. In an open chain there are (R+1-s) pairs at
    // separation s; in a closed chain (the (R+1)-th step lands back on the original) separations
    // s and (R+1-s) coincide, which the +1 on the first entry and the reduced loop bound encode.
    std::vector<CopyPair> out;
    bool closed = is_closed();
    int reps = static_cast<int>(repetitions());
    for (int k = 0; k < reps - static_cast<int>(closed); ++k) {
        int scale = reps - k;
        if (k == 0 && closed) {scale += 1;}
        out.push_back({0, k+1, scale});
    }
    return out;
}
