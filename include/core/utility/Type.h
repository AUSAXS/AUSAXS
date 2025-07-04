// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <string>
#include <typeinfo>

namespace ausaxs {
    /**
     * @brief Attempt to convert the internal type representation to a more human-readable format. 
     */
    std::string demangle(const char* name);

    /**
     * @brief Get the type of an object.
     *        This implementation is taken straight off stackoverflow: https://stackoverflow.com/questions/281818/unmangling-the-result-of-stdtype-infoname 
     */
    template <class T>
    std::string type(const T& t) {
        return demangle(typeid(t).name());
    }
}