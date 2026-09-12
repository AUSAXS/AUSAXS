// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <exception>
#include <string>

/**
 * @brief The exceptions available to the math subproject.
 */
namespace ausaxs::math::except {
    struct base : public std::exception {
        base(const char* msg);
        base(std::string msg);
        const char* what() const noexcept override;
        const std::string msg;
    };

    struct runtime_error : public base {using base::base;};
    struct invalid_argument : public base {using base::base;};
    struct out_of_range : public base {using base::base;};
    struct logic_error : public base {using base::base;};
    struct io_error : public base {using base::base;};
}
