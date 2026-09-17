// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <cstddef>

namespace ausaxs::table {
    struct DebyeTable {
        //! note: constexpr destructor cannot be defaulted due to GCC bug 93413
        constexpr virtual ~DebyeTable() noexcept = default;

        [[nodiscard]] virtual double lookup(int q_index, int d_index) const = 0;

        [[nodiscard]] virtual std::size_t size_q() const noexcept = 0;

        [[nodiscard]] virtual std::size_t size_d() const noexcept = 0;

        [[nodiscard]] virtual const double* begin(int q_index) const = 0;

        [[nodiscard]] virtual const double* end(int q_index) const = 0;
    };
}