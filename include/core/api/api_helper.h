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
        *status = 1;
    } catch (...) {
        ErrorMessage::last_error = "An unknown error occurred.";
        *status = 1;
    }
    if constexpr (!std::is_void_v<decltype(f())>) {return {};}
}