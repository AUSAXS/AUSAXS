// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <concepts>

namespace ausaxs {
    template<typename T>
        concept FormFactorType = requires(T t, double q) {
            {t.evaluate(q)};
    };
}