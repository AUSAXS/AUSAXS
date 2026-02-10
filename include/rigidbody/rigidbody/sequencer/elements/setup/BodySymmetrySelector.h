// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <cassert>

namespace ausaxs::rigidbody::sequencer::detail {
    struct BodySymmetrySelector {int body; int symmetry; int replica;};

    inline int to_index(int body, int symmetry=-1, int replica=0) {
        assert(symmetry >= -1 && symmetry < 100 && replica < 100);
        // symmetry=-1 indicates base body, so we add 1 to make it 0-indexed for storage
        return replica + (symmetry+1)*100 + body*10000;
    }

    inline BodySymmetrySelector from_index(int index) {
        assert(index < 1e6);
        int body = index / 10000;
        int symmetry = (index % 10000) / 100 - 1;  // subtract 1 to restore -1 for base bodies
        int replica = index % 100;
        return {.body=body, .symmetry=symmetry, .replica=replica};
    }
}