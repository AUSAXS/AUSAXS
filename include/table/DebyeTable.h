#pragma once

#include <constants/Axes.h>

namespace table {
    struct DebyeTable {
        // note: constexpr destructor cannot be defaulted due to GCC bug 93413
        constexpr virtual ~DebyeTable() noexcept {}

        [[nodiscard]] virtual constants::axes::d_type lookup(unsigned int q_index, unsigned int d_index) const = 0;

        [[nodiscard]] virtual unsigned int size_q() const noexcept = 0;

        [[nodiscard]] virtual unsigned int size_d() const noexcept = 0;

        [[nodiscard]] virtual const constants::axes::d_type* begin(unsigned int q_index) const = 0;

        [[nodiscard]] virtual const constants::axes::d_type* end(unsigned int q_index) const = 0;
    };
}