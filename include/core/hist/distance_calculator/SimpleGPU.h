#pragma once

// #define __WEBGPU__
#if defined(__WEBGPU__)
    #include <gpu/WebGPUSimple.h>
#elif defined(__ACPP__)
    #include <gpu/SimpleGPU.h>
#endif