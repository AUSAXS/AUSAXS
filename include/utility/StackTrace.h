#pragma once

#if defined(unix)
#include <backward.hpp>

namespace backward {
    backward::SignalHandling sh;
}
#endif