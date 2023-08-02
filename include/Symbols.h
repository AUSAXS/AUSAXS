#pragma once

#ifndef SAFE_MATH
    #define SAFE_MATH true
#endif

/**
 * @brief Print a message if DEBUG is enabled.
 */
#ifdef DEBUG
    #define debug_print(msg) \
        do { std::cout << msg << std::endl; } while (0)
#else
    #define debug_print(msg) \
        do { } while (0)
#endif
