#pragma once

#if defined(__ACPP__)
    #include <gpu/SimpleGPU.h>
#elif defined(__WEBGPU__)
    #include <gpu/WebGPUSimple.h>
#endif