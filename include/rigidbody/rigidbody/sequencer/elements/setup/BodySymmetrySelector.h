// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <cassert>

namespace ausaxs::rigidbody::sequencer::detail {
    struct BodySymmetrySelector {int body; int symmetry; int replica;};

    inline int to_index(int body, int symmetry=0, int replica=0) {
        assert(symmetry < 100 && replica < 100);
        return replica + symmetry*100 + body*10000;
    }

    inline BodySymmetrySelector from_index(int index) {
        assert(index < 1e6);
        int body = index / 10000;
        int symmetry = (index % 10000) / 100;
        int replica = index % 100;
        return {.body=body, .symmetry=symmetry, .replica=replica};
    }
}