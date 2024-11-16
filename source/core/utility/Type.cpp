/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <utility/Type.h>
#ifdef __GNUG__
#include <cstdlib>
#include <memory>
#include <cxxabi.h>

using namespace ausaxs;

std::string demangle(const char* name) {
    int status = -4; // some arbitrary value to eliminate the compiler warning
    std::unique_ptr<char, void(*)(void*)> res {abi::__cxa_demangle(name, NULL, NULL, &status), std::free};
    return (status==0) ? res.get() : name ;
}

#else

// does nothing if not g++
std::string demangle(const char* name) {return name;}

#endif