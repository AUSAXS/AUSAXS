#pragma once

#ifdef MSVC
    #define API __declspec(dllexport)
#else
    #define API
#endif

extern "C" API void evaluate_sans_debye(double* _q, double* _x, double* _y, double* _z, double* _w, int _nq, int _nc, int _return_status, double* _return_Iq);