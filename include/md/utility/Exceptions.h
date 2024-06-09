#pragma once

#include <string>

namespace except {
    struct base : public std::exception {
        base(const char* msg) : msg(msg) {}
        base(const std::string msg) : msg(msg) {}
        const char* what() const noexcept {return msg.data();}
        const std::string msg;
    };

    struct invalid_order : public base {using base::base;};
    struct unexpected : public base {using base::base;};
    struct duplicate_option : public base {using base::base;};
    struct missing_option : public base {using base::base;};
    struct invalid_argument : public base {using base::base;};
    struct invalid_format : public base {using base::base;};
    struct io_error : public base {using base::base;};
    struct unknown_type : public base {using base::base;};
}