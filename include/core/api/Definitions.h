#pragma once

#ifdef WIN32
    #define API __declspec(dllexport)
#else
    #define API
#endif