// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

namespace ausaxs {
    /**
     * @brief An observer-pointer is just a wrapper around a raw pointer to denote that ownership is not transferred.
     *        This is conceptually equivalent to std::experimental::observer_ptr.
     */
    template<typename T> using observer_ptr = T*;
}