#include <api/ErrorMessage.h>

extern "C" API void get_last_error_msg(
    char** buffer, int* buffer_size, 
    int* status
) {
    *status = 1;
    *buffer_size = static_cast<int>(ErrorMessage::last_error.size())+1;
    *buffer = const_cast<char*>(ErrorMessage::last_error.c_str());
    *status = 0;
}