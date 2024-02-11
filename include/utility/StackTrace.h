#pragma once

#if defined(unix) && defined(DEBUG)
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