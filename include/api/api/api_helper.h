// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <string>
#include <type_traits>

#ifdef WIN32
    #define API __declspec(dllexport)
#else
    #define API
#endif

struct ErrorMessage {
    inline static std::string last_error;
};

extern "C" API void get_last_error_msg(
    char** msg,
    int* status
);

/**
 * @brief Write an API-level exception to the error stream.
 */
void report_api_exception(const char* what) noexcept;

template<typename Func>
auto execute_with_catch(Func&& f, int* status) -> decltype(f()) {
    try {
        *status = 1;
        if constexpr (std::is_void_v<decltype(f())>) {
            f();
            *status = 0;
            return;
        } else {
            auto v = f();
            *status = 0;
            return v;
        }
    } catch (const std::exception& e) {
        ErrorMessage::last_error = e.what();
        report_api_exception(ErrorMessage::last_error.c_str());
        *status = 1;
    } catch (...) {
        // not derived from std::exception, so there is nothing to report but its existence
        ErrorMessage::last_error = "An unknown error occurred.";
        report_api_exception(ErrorMessage::last_error.c_str());
        *status = 1;
    }
    if constexpr (!std::is_void_v<decltype(f())>) {return {};}
}