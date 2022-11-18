#pragma once

#ifndef SAFE_MATH
    #define SAFE_MATH true
#endif

#ifndef DEBUG 
    #define DEBUG true
#endif

/**
 * @brief Print a message if DEBUG is enabled.
 */
#define debug_print(msg) \
    do { if (DEBUG) std::cout << msg << std::endl; } while (0)