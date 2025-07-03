// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#if defined(DEBUG)
    #define BACKWARD_HAS_DWARF 1
    #include <backward.hpp>

    namespace backward {
        static inline backward::SignalHandling sh;
    }

    inline void print_trace() {
        backward::StackTrace st;
        st.load_here(32);
        backward::Printer p;
        p.print(st);
    }
#endif