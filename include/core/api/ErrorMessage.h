#pragma once

#include <api/Definitions.h>

#include <string>

struct ErrorMessage {
    inline static std::string last_error;
};

extern "C" API void get_last_error_msg(
    char** msg,
    int* status
);
