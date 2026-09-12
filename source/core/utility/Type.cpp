// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <utility/Type.h>

#ifdef __GNUG__
    #include <cstdlib>
    #include <cxxabi.h>
    #include <memory>

    std::string ausaxs::demangle(const char* name) {
        int status = -4; // some arbitrary value to eliminate the compiler warning
        std::unique_ptr<char, void(*)(void*)> res {abi::__cxa_demangle(name, nullptr, nullptr, &status), std::free};
        return (status==0) ? res.get() : name ;
    }
#else
    // does nothing if not g++
    std::string ausaxs::demangle(const char* name) {return name;}
#endif