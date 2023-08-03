#pragma once

#if defined(unix) && defined(DEBUG)
    #include <backward.hpp>

    namespace backward {
        backward::SignalHandling sh;
    }
#endif