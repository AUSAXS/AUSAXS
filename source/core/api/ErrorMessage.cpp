#include <api/ErrorMessage.h>

extern "C" API void get_last_error_msg(
    char** buffer,
    int* status
) {
    *status = 1;
    *buffer = const_cast<char*>(ErrorMessage::last_error.c_str());
    *status = 0;
}